import pandas as pd
import os
from pathlib import Path
import argparse
import toml
import sys
import subprocess

try:
    parser = argparse.ArgumentParser(
        prog="Rename Bam files after dorado demultiplexing",
        description="Change long names for bam files to a easier format",
    )

    parser.add_argument(
        "--config", type=str, required=True, help="Project config file, including path."
    )
    parser.add_argument("--sample", type=str, required=True, help="Sample.")

    args = parser.parse_args()

    # Loading TOML config
    with open(args.config, "r") as f:
        toml_config = toml.load(f)

    flowcells = toml_config["general"]["fc_dir_names"]

    username = os.environ.get("USER")
    tmpdir = os.environ.get("SLURM_TMPDIR")
    tmp = Path(tmpdir)
    dir_proj = toml_config["general"]["project_path"]
    name = dir_proj.rstrip("/").split("/")[-2].split("_", 1)[1]
    output = f"/lustre10/scratch/{username}/{name}/alignments"
    output_path = Path(output)

    s = args.sample

    all_inputs = []

    # RENAME
    for fc in flowcells:
        print(f"\nRunning: rename for flowcell {fc}")
        code = fc.split("_")[-1]

        inputs = f"/lustre10/scratch/{username}/{name}/{fc}/"
        inputs = Path(inputs)
        all_inputs.append(inputs)

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

                        if new_name.startswith(s):
                            new_path = inputs / new_name
                            file.rename(new_path)
                            print(f"Renamed {file.name} -> {new_name}")

    # SAMTOOLS
    print(f"\nRunning Samtools for sample {s}")
    output_file = tmp / f"{s}.bam"

    bam_files = []
    for p in all_inputs:
        bam_files.extend(p.glob(f"{s}_*.bam"))

    bam_files_str = [str(p) for p in bam_files]

    # merge
    print("Merge")
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

    # transfer bam
    transfer1 = ["cp", str(output_file), output_path / f"{s}.bam"]
    subprocess.run(transfer1, check=True)
    rm1 = ["rm", str(output_file)]
    subprocess.run(rm1, check=True)

    # sort
    print("Sort")
    bam = output_path / f"{s}.bam"
    sorted_bam = tmp / f"{s}_sorted.bam"
    cmd2 = [
        "samtools",
        "sort",
        "--threads",
        "3",
        "-m",
        "4G",
        "-o",
        str(sorted_bam),
        str(bam),
    ]
    subprocess.run(cmd2, check=True)

    # transfer 2
    transfer2 = ["cp", str(sorted_bam), str(output_path / f"{s}_sorted.bam")]
    subprocess.run(transfer2, check=True)

    # index
    print("Index")
    cmd3 = ["samtools", "index", "--threads", "3", str(output_path / f"{s}_sorted.bam")]
    subprocess.run(cmd3, check=True)

    print(f"Done Samtools for {s}")

except Exception as e:
    print(f"Error: {e}")
    sys.exit(1)
