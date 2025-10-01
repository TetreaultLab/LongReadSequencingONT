import pandas as pd
from pathlib import Path
from collections import defaultdict
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

    args = parser.parse_args()

    # Loading TOML config
    with open(args.config, "r") as f:
        toml_config = toml.load(f)

    flowcells = ast.literal_eval(args.flowcells)

    all_inputs = []
    # Loop over flowcells to rename
    for fc in flowcells:
        print("\nRunning: rename for flowcell ", fc)
        code = fc.split("_")[-1]

        inputs = toml_config["general"]["project_path"] + "/" + fc + "/alignments"
        inputs = Path(inputs)
        all_inputs.append(inputs)

        output = toml_config["general"]["project_path"] + "/alignments"
        output_path = Path(output)

        # Load the CSV file
        df = pd.read_csv(
            toml_config["general"]["project_path"] + "/scripts/" + fc + ".csv", header=0
        )
        df["code"] = df["flow_cell_id"].str.split("_").str[-1]

        # Create a mapping from barcode -> alias
        barcode_to_alias = dict(zip(zip(df["barcode"], df["alias"]), df["code"]))

        # Process each file in the directory
        for file in inputs.iterdir():
            if file.is_file() and file.name.startswith(code):
                for (barcode, alias), code in barcode_to_alias.items():
                    if barcode in file.name:
                        if file.name.endswith(".bam"):
                            new_name = f"{alias}_{barcode}_{code}.bam"
                        else:
                            continue

                        new_path = inputs / new_name
                        file.rename(new_path)
                        print(f"Renamed {file.name} -> {new_name}")

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

        print(f"Done {s}")

    # Run all samples in parallel
    samples = toml_config["general"]["samples"]
    n = len(samples)
    with Pool(processes=n) as pool:
        pool.map(process_sample, samples)


except Exception as e:
    print(f"Error: {e}")
    sys.exit(1)
