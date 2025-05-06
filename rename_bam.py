import pandas as pd
from pathlib import Path
import argparse

parser = argparse.ArgumentParser(description="Rename BAM and BAI files based on barcode->alias mapping.")
parser.add_argument("directory", type=str, help="Path to the directory containing the BAM/BAI files.")

# Load the CSV file
df = pd.read_csv(args.directory + "/samples.csv")

# Create a mapping from barcode -> alias
barcode_to_alias = dict(zip(df["barcode"], df["alias"]))

# Directory containing the files
data_dir = Path(args.directory + "/alignments")

# Process each file in the directory
for file in data_dir.iterdir():
    if file.is_file():
        for barcode, alias in barcode_to_alias.items():
            if barcode in file.name:
                if file.name.endswith(".bam.bai"):
                    new_name = f"{alias}_{barcode}.bam.bai"
                else:
                    new_name = f"{alias}_{barcode}.bam"
                    
                new_path = data_dir / new_name
                file.rename(new_path)
                print(f"Renamed {file.name} -> {new_name}")
