import pandas as pd
from pathlib import Path
import argparse
import toml
import sys
import os

try:
    parser = argparse.ArgumentParser(
        prog="Rename Bam files after dorado demultiplexing",
        description="Change long names for bam files to a easier format",
    )

    parser.add_argument(
        "--config", type=str, required=True, help="Project config file, including path."
    )
    parser.add_argument("--flowcell", type=str, required=True, help="Flowcell.")

    args = parser.parse_args()

    fc = args.flowcell

    print(f"\nRunning: rename for flowcell {fc}")

    # Loading TOML config
    with open(args.config, "r") as f:
        toml_config = toml.load(f)

    code = fc.split("_")[-1]

    slurmtmp = os.environ.get("SLURM_TMPDIR")
    inputs = f"{slurmtmp}/{fc}"
    inputs = Path(inputs)

    # Load the CSV file
    df = pd.read_csv(
        toml_config["general"]["project_path"] + "/scripts/" + fc + ".csv", header=0
    )
    df["code"] = df["flow_cell_id"].str.split("_").str[-1]

    # Create a mapping from barcode -> alias
    barcode_to_alias = dict(zip(zip(df["barcode"], df["alias"]), df["code"]))

    # Process each file in the directory
    for file in inputs.iterdir():
        print(file)
        if file.is_file() and file.name.startswith(code):
            for (barcode, alias), code in barcode_to_alias.items():
                if barcode in file.name:
                    if file.name.endswith(".bam"):
                        print(file.name)
                        new_name = f"{alias}_{barcode}_{code}.bam"
                    else:
                        print(file.name)
                        continue

                    new_path = inputs / new_name
                    # file.rename(new_path)
                    print(f"Renamed {file.name} -> {new_name}")

    print(f"\nDone: rename for flowcell {fc}")

except Exception as e:
    print(f"Error: {e}")
    sys.exit(1)
