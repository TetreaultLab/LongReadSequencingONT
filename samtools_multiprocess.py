import pandas as pd
import os
from pathlib import Path
import argparse
import toml
import subprocess
import ast
import sys
from multiprocessing import Pool

try:
    parser = argparse.ArgumentParser(
        prog="Rename Bam files after dorado demultiplexing",
        description="Change long names for bam files to a easier format",
    )

    parser.add_argument(
        "--config", type=str, required=True, help="Project config file, including path."
    )
    parser.add_argument(
        "--flowcells", type=str, required=True, help="Flowcells to run."
    )
    parser.add_argument("--sample", type=str, required=True, help="Sample.")

    args = parser.parse_args()

    # Loading TOML config
    with open(args.config, "r") as f:
        toml_config = toml.load(f)

    flowcells = ast.literal_eval(args.flowcells)

    username = os.environ.get("USER")
    dir_proj = toml_config["general"]["project_path"]
    name = dir_proj.rstrip("/").split("/")[-2].split("_", 1)[1]

    output = f"/lustre10/scratch/{username}/{name}/alignments"
    output_path = Path(output)

    all_inputs = []
    # Loop over flowcells to get inputs
    for fc in flowcells:
        inputs = f"/lustre10/scratch/{username}/{name}/{fc}/"
        inputs = Path(inputs)
        all_inputs.append(inputs)

    s = args.sample

    def process_sample(s):
        print(f"\nRunning Samtools for sample {s}")
        output_file = output_path / f"{s}.bam"

        bam_files = []
        for p in all_inputs:
            bam_files.extend(p.glob(f"{s}_*.bam"))

        bam_files_str = [str(p) for p in bam_files]

        # merge
        cmd = [
            "samtools",
            "merge",
            "-c",
            "-f",
            "--threads",
            "3",
            "-o",
            str(output_file),
        ] + bam_files_str
        subprocess.run(cmd, check=True)

        # sort
        sorted_bam = output_path / f"{s}_sorted.bam"
        cmd2 = [
            "samtools",
            "sort",
            "--threads",
            "3",
            "-m",
            "4G",
            "-o",
            str(sorted_bam),
            str(output_file),
        ]
        subprocess.run(cmd2, check=True)

        # index
        cmd3 = ["samtools", "index", "--threads", "3", str(sorted_bam)]
        subprocess.run(cmd3, check=True)

        print(f"Done Samtools for {s}")

except Exception as e:
    print(f"Error: {e}")
    sys.exit(1)
