import argparse
import toml
import os
from glob import glob

parser = argparse.ArgumentParser(
    prog="PipelineLongReads",
    description="LongReadSequencing pipeline.",
    )

parser.add_argument(
    "config", type=str, help="Project config file, including path."
    )

args = parser.parse_args()

# Read initial config file
with open(args.config, "r") as f:
    toml_config = toml.load(f)
    print("~~~ INITIAL CONFIG ~~~")
    print(toml_config)

## Which analysis are possible for each seq-type. Make a if to fill the analysis field and check if some error are made in analysis was not empty
# if analysis list is empty, add all available analysis
if toml_config["general"]["analysis"] == []:
    toml_config["general"]["analysis"] = ["methylation", "splicing", "polya", "snp", "sv", "repeats"]

# set kit to variable
kit = toml_config['general']['kit']

# check seq_type
seq_type = toml_config["general"]["seq_type"]
if (toml_config["general"]["seq_type"] == "Meth") or ("methylation" in toml_config["general"]["analysis"]):
    methylation_status = True
else:
    methylation_status = False

# if path has a "/" at the end, remove it 
path = toml_config["general"]["project_path"]
if path.endswith('/'):
    new_path = path[0:-1]
    toml_config["general"]["project_path"] = new_path
else:
    new_path = path

# Making directory structure
directories = ["main_reports", "reads", "scripts", "alignments", "results", "qc"]
for d in directories:
    if not os.path.exists(new_path + "/" + d):
        os.makedirs(new_path + "/" + d)

# Move main reports files to corresponding directory
reports = ["barcode_alignment", "final_summary", "pore_activity", "report", "sample_sheet", "sequencing_summary", "throughput"]
for r in reports:
    f = glob(os.path.join(new_path, r + "*"))
    if f != []: # if the file exist in the main directory
        for i in range(0, len(f)): # if one or more file starts with the report name
            f2 = f[i]
            name = f2.split("/")[-1]
            os.rename(f2, os.path.join(new_path, "main_reports", name)) # move to main_reports

# Add Dorado options
## general options
toml_config["dorado"] = {}
toml_config['dorado']['min_q_score'] = 10
toml_config['dorado']['sample_sheet'] = "samples.csv"
toml_config['dorado']['barcode_both_ends'] = True
toml_config['dorado']['trim'] = "all"

## parameters depending on sequencing type
if seq_type == "WGS":
    toml_config['dorado']['mm2_opts'] = "-ax lr:hq"
    if kit not in ["SQK-RBK114-24", "SQK-NBD114.24", "SQK-LSK114"]:
        raise Exception("Error: Wrong Kit for WGS. Options are SQK-RBK114-24, SQK-NBD114.24, SQK-LSK114")

if seq_type == "RNA":
    toml_config['dorado']['mm2_opts'] = "-ax splice:hq -uf"
    if kit not in ["SQK-PCB114-24"]:
        raise Exception("Error: Wrong Kit for Whole Transcriptome. Options are SQK-PCB114-24")

if methylation_status:
    toml_config['dorado']['modified_bases'] = "5mCG_5hmCG"
    toml_config['dorado']['modified_bases_threshold'] = 0.05

if seq_type == "Targeted":
    toml_config['dorado']['mm2_opts'] = "-ax splice --junc-bed anno.bed12"
    if kit not in ["SQK-NBD114-24"]:
        raise Exception("Error: Wrong Kit for Targeted Sequencing. Options are SQK-NBD114-24")


# Add next tool options


# Save new config file
print("~~~ FINAL CONFIG ~~~")
print(toml_config)
with open(toml_config["general"]["project_path"] + '/config_final.toml', 'w') as f:
    toml.dump(toml_config, f)

