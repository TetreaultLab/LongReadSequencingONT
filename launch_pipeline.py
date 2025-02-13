import argparse
import toml
from datetime import datetime
import subprocess
import sys
import pandas as pd
import numpy as np

def main():
    start = get_time()
    
    parser = argparse.ArgumentParser(
        prog="PipelineLong",
        description="LongReadSequencing pipeline.",
    )

    parser.add_argument(
        "--config", type=str, required=True, help="Project config file, including path."
    )

    args = parser.parse_args()
    
    # Loading TOML config
    with open(args.config, "r") as f:
        toml_config = toml.load(f)

    # Start of pipeline
    start_str = ">>> LRS pipeline starting at {}.".format(
        start
    )
    print(
        "=" * len(start_str) + "\n" + start_str + "\n" + "=" * len(start_str),
        file=sys.stdout,
    )
    
    output = toml_config["general"]["project_path"]

    # Open file for steps done
    steps = open(output + "/steps_done.txt", "a")
    steps.write("\nLoading ENV\n")
    steps.close()
    
    # Create list of steps already done from the file steps_done.txt
    done = []
    with open(output + "/steps_done.txt", "r") as f:
        for line in f:
            done.append(line.strip())


    # Get tools and versions
    function_queue = []

    genome = get_reference(toml_config["general"]["reference"], "")["fasta"]
    print(f"\t>>> Reference genome version: {genome}")
    
    # Setting up list of steps
    # Base calling
    if "dorado" not in done:
        print("\t>>> Base calling: Dorado (?)")
        function_queue.append(dorado)

    # Other tools ...


    # Calling each steps
    for func in function_queue:
        try:
            func(toml_config)
        except:
            exit(1)

    end = get_time()
    total_time = end - start
    end_str = ">>> LRS-seq pipeline completed in {}.".format(
        total_time
    )
    print("=" * len(end_str) + "\n" + end_str + "\n" + "=" * len(end_str))


def get_time():
    now = datetime.now()
    return now


def title(message):
    print(f"\n>>> Running {message} [{get_time()}]")


### TO CHANGE
def get_reference(ref, tool):
    path = "/lustre03/project/6019267/shared/tools/PIPELINES/References/"
    reference: {}  # type: ignore
    match ref:
        case "grch38":
            reference = {
                "fasta": path + "Homo_sapiens.GRCh37.87.dna.primary_assembly.fa",
                "index": path + "index_" + tool + "/" + ref,
                "gtf": path + "Homo_sapiens.GRCh37.87.gtf",
                "gff3": path + "Homo_sapiens.GRCh37.87.gff3.gz",
            }

    return reference


def dorado(toml_config):
    



if __name__ == "__main__":
    main()
