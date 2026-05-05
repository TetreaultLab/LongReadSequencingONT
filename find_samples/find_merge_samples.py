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

    print(f"########## Project name : {project} ##########")
    print(print(*samples, sep="\n"))

    output = (
        f"/lustre09/project/6019267/shared/projects/Nanopore_Dock/Combined/{project}"
    )
    dir_path = Path(output)

    # check if project name is unique and create output directory
    Path(dir_path).mkdir(exist_ok=True)

    # for each sample run the find function to create a file with bam paths.
    for s in samples:
        find(s)


def find(sample):
    print(f">>> Looking for {sample}")
    for path in Path("/lustre09/project/6019267/shared/projects/Nanopore_Dock").rglob(
        f"{sample}.bam"
    ):
        print(path)


def merge():
    print("Merge")


# Launches main function
if __name__ == "__main__":
    main()
