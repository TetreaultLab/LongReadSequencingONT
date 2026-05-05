#!/usr/bin/env python3

import argparse
from pathlib import Path


def main():
    print("Executing main\n")

    parser = argparse.ArgumentParser(
        prog="Find samples in Nanopore_Dock",
        description="Pipeline to compare LRS samples.",
    )

    parser.add_argument(
        "--project", type=str, required=True, help="Project name. Must be unique"
    )

    group = parser.add_mutually_exclusive_group(required=True)

    group.add_argument("--rna", action="store_true", help="RNA-seq analysis")
    group.add_argument("--genome", action="store_true", help="WGS analysis")

    parser.add_argument(
        "--samples",
        nargs="+",
        type=str,
        required=True,
        help="Samples to include in project. Space delimited, e.g.: --samples S1 S2 S3",
    )

    args = parser.parse_args()

    project = args.project
    samples = args.samples
    if args.rna:
        seq = ["RNA", "rna", "Transcriptome", "transcriptome"]
    if args.genome:
        seq = ["WGS", "wgs"]

    print(f"\n########## Project name : {project} ##########")
    print("Looking for sequencing type :")
    print(*seq, sep="/")
    print("For samples :")
    print(*samples, sep="\n")

    output = (
        f"/lustre09/project/6019267/shared/projects/Nanopore_Dock/Combined/{project}"
    )
    dir_path = Path(output)

    # check if project name is unique and create output directory
    Path(dir_path).mkdir(exist_ok=False)

    # for each sample run the find function to create a file with bam paths.
    for s in samples:
        find(s)


def find(sample, seq):
    print(f"\n>>> Looking for {sample}")
    nanopore = Path("/lustre09/project/6019267/shared/projects/Nanopore_Dock/")
    # Look only in directories starting with a "2" and not containing "Unified" and not in other directories as they could already be combined
    for subdir in nanopore.iterdir():
        if (
            subdir.is_dir()
            and subdir.name.startswith("2")
            and "Unified" not in subdir.name
        ):
            # Recursively search within each matching directory
            paths = list(subdir.rglob(f"{sample}_sorted.bam"))
    print(paths)

    # Save paths to file
    with open(f"{str(nanopore)}/{sample}_bam_paths.txt", "w") as f:
        for p in paths:
            f.write(f"{str(p)}\n")


def merge():
    print("Merge")


def create_config():
    print(
        "Create confing for merged samples and launch pipeline for downstream analyses"
    )


# Launches main function
if __name__ == "__main__":
    main()
