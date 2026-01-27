import argparse
import toml
import sys
import time
import subprocess
import pandas as pd
import numpy as np
import os
from glob import glob

# Global variables
TOOL_PATH = "/lustre09/project/6019267/shared/tools/"
OUTPUT = "/lustre09/project/6019267/shared/projects/GenomeCanada_CPHI_NRGI/results/LRS/comparisons/"


def main():
    parser = argparse.ArgumentParser(
        prog="PipelineLong",
        description="LongReadSequencing pipeline.",
    )

    parser.add_argument("config", type=str, help="Project config file.")

    args = parser.parse_args()

    # Loading initial TOML config
    with open(args.config, "r") as f:
        toml_config = toml.load(f)

    # check if comparison directory alredy exists
    folder_name = OUTPUT + toml_config["general"]["comparison_name"]

    if os.path.exists(folder_name):
        print(
            f"The path '{folder_name}' exists. \nExiting. \nPlease choose a new 'comparison_name' and relaunch."
        )
        sys.exit(1)

    else:
        print(f"Making the folder '{folder_name}'.")
        try:
            os.makedirs(folder_name)
            print(f"Success! The folder '{folder_name}' was created.")

        except OSError as e:
            print(f"Error creating directory '{folder_name}': {e}")
            sys.exit(1)

    # Get functions to run
    function_queue = []

    if "degs" in toml_config["general"]["analyses"]:
        function_queue.append(flair_diffexp)

    if "as" in toml_config["general"]["analyses"]:
        function_queue.append(flair_diffsplice)

    if "apa" in toml_config["general"]["analyses"]:
        function_queue.append(apa_tool)

    if "dmrs" in toml_config["general"]["analyses"]:
        function_queue.append(modkit)

    # Calling the functions
    for func in function_queue:
        func(toml_config)


def flair_diffexp(toml_config):
    print()


def flair_diffsplice(toml_config):
    print()


def apa_tool(toml_config):
    print()


def modkit(toml_config):
    print()


# Launches main function
if __name__ == "__main__":
    main()
