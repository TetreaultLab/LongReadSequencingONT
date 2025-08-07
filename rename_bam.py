import pandas as pd
from pathlib import Path
import toml

# Loading TOML config
with open(args.config, "r") as f:
    toml_config_initial = toml.load(f)

output = toml_config_initial["general"]["project_path"] + "/alignments"
fcs = toml_config_initial["general"]["fc_dir_names"]

# Loop over flowcells to rename
for fc in fcs :
    code = fc.split('_')[-1]
    print(code)

    # Load the CSV file
    df = pd.read_csv(toml_config_initial["general"]["project_path"] + "/scripts/" + fc + ".csv", header=0)
    df["code"] = df["flow_cell_id"].str.split('_').str[-1]
    print(df)

    # Create a mapping from barcode -> alias
    barcode_to_alias = dict(zip(df["barcode"], df["alias"], df["code"]))

    # Process each file in the directory
    for file in output.iterdir():
        if file.is_file() and file.name.startswith(code):
            for barcode, alias, code in barcode_to_alias.items():
                if barcode in file.name:
                    if file.name.endswith(".bam.bai"):
                        new_name = f"{alias}_{barcode}_{code}.bam.bai"
                    else:
                        new_name = f"{alias}_{barcode}_{code}.bam"
                        
                    new_path = output / new_name
                    #file.rename(new_path) # uncomment to actually change the names
                    print(f"Renamed {file.name} -> {new_name}")

# Merge bams
# TO-DO
