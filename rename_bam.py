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
    toml_config = toml.load(f)

fcs = toml_config["general"]["fc_dir_names"]

all_inputs = []
# Loop over flowcells to rename
for fc in fcs :
    print("\nRunning: rename for flowcell ", fc)
    code = fc.split('_')[-1]

    inputs = toml_config["general"]["project_path"] + "/" + fc + "/alignments"
    inputs = Path(inputs)
    all_inputs.append(inputs)

    output = toml_config["general"]["project_path"] + "/alignments"
    output_path = Path(output)

    # Load the CSV file
    df = pd.read_csv(toml_config["general"]["project_path"] + "/scripts/" + fc + ".csv", header=0)
    df["code"] = df["flow_cell_id"].str.split('_').str[-1]

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


# Merge bams
samples = toml_config["general"]["samples"]

# Loop over samples
for s in samples:
    print("\nRunning: Samtools for sample ", s)
    output_file = output_path / f"{s}.bam"

    bam_files = []
    for p in all_inputs:
        bam_files.extend(p.glob(f"{s}_*.bam"))
    print(bam_files)
    
    bam_files_str = ", ".join(str(f) for f in bam_files)
    print(bam_files_str)

    # merge
    print("--> Merge")
    cmd = ["samtools", "merge", "-f", "-@", "8", "-o", str(output_file)] + bam_files
    print(cmd)
    subprocess.run(cmd, check=True)

    # sort
    print("--> Sort")
    cmd2 = ["samtools", "sort", "-@", "8", "-o", output + "/" + s + "_sorted.bam", output + "/" + s + ".bam"]
    subprocess.run(cmd2, check=True)
    
    # index
    print("--> Index")
    cmd3 = ["samtools", "index", "-@", "8", "-o", output + "/" + s + "_sorted.bam.bai", output + "/" + s + "_sorted.bam"]
    subprocess.run(cmd3, check=True)

path = toml_config["general"]["project_path"]
rm_prefix = path.replace('/lustre09/project/6019267/shared/projects/Nanopore_Dock/', '')
#rm_prefix = path.replace('/lustre09/project/6019267/shared/tools/main_pipelines/long-read/', '')
path_list = rm_prefix.split("/")
project_name_date = path_list[0].split("_", 1)
project_name = project_name_date[1]
print("Done ", project_name, " !")

