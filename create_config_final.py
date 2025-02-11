import argparse
import toml

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
    print(kit)
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

