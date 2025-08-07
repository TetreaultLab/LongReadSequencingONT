import pandas as pd
from pathlib import Path
from collections import defaultdict
import argparse
import toml
import subprocess

parser = argparse.ArgumentParser(
        prog="Rename Bam files after dorado demultiplexing",
        description="Change long names for bam and .bam.bai files to a easier format",
    )

parser.add_argument(
    "config", type=str, help="Project config file, including path."
    )

args = parser.parse_args()

# Loading TOML config
with open(args.config, "r") as f:
    toml_config_initial = toml.load(f)

output = toml_config_initial["general"]["project_path"] + "/alignments"
output = Path(output)
fcs = toml_config_initial["general"]["fc_dir_names"]

# Loop over flowcells to rename
for fc in fcs :
    code = fc.split('_')[-1]

    # Load the CSV file
    df = pd.read_csv(toml_config_initial["general"]["project_path"] + "/scripts/" + fc + ".csv", header=0)
    df["code"] = df["flow_cell_id"].str.split('_').str[-1]

    # Create a mapping from barcode -> alias
    barcode_to_alias = dict(zip(zip(df["barcode"], df["alias"]), df["code"]))

    # Process each file in the directory
    for file in output.iterdir():
        if file.is_file() and file.name.startswith(code):
            for (barcode, alias), code in barcode_to_alias.items():
                if barcode in file.name:
                    if file.name.endswith(".bam"):
                        new_name = f"{alias}_{barcode}_{code}.bam"
                    else:
                        continue
                        
                    new_path = output / new_name
                    file.rename(new_path)
                    print(f"Renamed {file.name} -> {new_name}")


# Merge bams
samples = toml_config_initial["general"]["samples"]

# Loop over samples
for s in samples:
    print("Running: samtools merge, sort and index for sample ", s)
    bam_files = list(output.glob(f"{s}*.bam"))
    output_file = output / f"{s}.bam"

    # merge
    cmd = ["samtools", "merge", "-o", str(output_file)] + [str(f) for f in bam_files]
    subprocess.run(cmd, check=True)

    # sort
    cmd2 = ["samtools", "sort", "-o", s + "_sorted.bam", s + ".bam"]
    subprocess.run(cmd2, check=True)
    
    # index
    cmd3 = ["samtools", "index", "-o", s + "_sorted.bam.bai", s + "_sorted.bam"]
    subprocess.run(cmd3, check=True)

# Clean-up : TO-DO

