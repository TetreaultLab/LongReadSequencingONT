import pandas as pd
from pathlib import Path
from collections import defaultdict
import argparse
import toml
import subprocess

parser = argparse.ArgumentParser(
        prog="Rename Bam files after dorado demultiplexing",
        description="Change long names for bam files to a easier format",
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
    print("\nRunning: rename for flowcell ", fc)
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
    print("Running: Samtools for sample ", s)
    bam_files = list(output.glob(f"{s}_*.bam"))
    output_file = output / f"{s}.bam"

    # merge
    print("--> Merge")
    cmd = ["samtools", "merge", "-@", "8", "-o", str(output_file)] + [str(f) for f in bam_files]
    subprocess.run(cmd, check=True)

    # sort
    print("--> Sort")
    cmd2 = ["samtools", "sort", "-@", "8", "-o", s + "_sorted.bam", s + ".bam"]
    subprocess.run(cmd2, check=True)
    
    # index
    print("--> Index")
    cmd3 = ["samtools", "index", "-@", "8", "-o", s + "_sorted.bam.bai", s + "_sorted.bam"]
    subprocess.run(cmd3, check=True)

path = toml_config["general"]["project_path"]
rm_prefix = path.replace('/lustre09/project/6019267/shared/tools/main_pipelines/long-read/', '')
path_list = rm_prefix.split("/")
project_name_date = path_list[0].split("_", 1)
project_name = project_name_date[1]
print("Done ", project_name, " !")

# Clean-up : TO-DO

