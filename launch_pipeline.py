import argparse
import toml
from datetime import datetime
import subprocess
import sys
import pandas as pd
import numpy as np
import traceback

def main():
    start = get_time()
    
    parser = argparse.ArgumentParser(
        prog="PipelineLong",
        description="LongReadSequencing pipeline.",
    )

    parser.add_argument(
        "config", type=str, help="Project config file, including path."
    )

    args = parser.parse_args()
    
    # Loading TOML config
    with open(args.config, "r") as f:
        toml_config = toml.load(f)
    
    print(toml_config)

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
    print(f">>> Reference genome version: {genome}")
    
    # Setting up list of steps
    # Base calling
    if "dorado" not in done:
        print(">>> Base calling: Dorado (?)")
        function_queue.append(dorado)

    # Other tools ...


    # Calling each steps
    for func in function_queue:
        func(toml_config)
            

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


def saving(toml_config, tool):
    print(f"\n>>> Finished running {tool} successfully [{get_time()}]")
    with open(toml_config["general"]["project_path"] + "/steps_done.txt", "a") as steps:
        steps.write(tool)
        steps.write("\n")


### TO CHANGE
def get_reference(ref, tool):
    path = "/home/shared/tools/reference_files/"
    reference: {}  # type: ignore
    match ref:
        case "grch38":
            reference = {
                "fasta": path + "hg38/gencode.v38.p13.genome.fa",
                "index": path + "index_" + tool + "/" + ref,
                "gtf": path + "Homo_sapiens.GRCh37.87.gtf",
                "gff3": path + "Homo_sapiens.GRCh37.87.gff3.gz",
            }

    return reference


def dorado(toml_config):
    title("dorado")

    output = toml_config["general"]["project_path"]
    genome = get_reference(toml_config["general"]["reference"], "dorado")["fasta"]
    reads = "/home/shared/data/2025-01-15_FXN-Batch4/FRDA14_21-UTMAB-06_2/20250115_2147_P2S-02441-B_PBA20836_7fce705b/pod5/PBA20836_7fce705b_9c89ba7b_66.pod5"
    
    command = ["dorado", "basecaller", "--verbose", "--device", "cuda:0", "--min-qscore", str(toml_config["dorado"]["min_q_score"]), "-o", output, "--reference", genome, "--sample-sheet", output + "/" + toml_config["dorado"]["sample_sheet"], "--trim", toml_config["dorado"]["trim"], "--kit-name", toml_config["general"]["kit"], "--mm2-opts", toml_config["dorado"]["mm2_opts"]]
    
    if toml_config["dorado"]["barcode_both_ends"] in ["true", "True", "yes", "Yes"]:
        command.extend(["--barcode-both-ends"])

    if "methylation" in toml_config["general"]["analysis"]:
        command.extend(["--modified-bases", toml_config["dorado"]["modified_bases"], "--modified-bases-threshold", str(toml_config["dorado"]["modified_bases_threshold"])])

    if "polya" in toml_config["general"]["analysis"]:
        command.extend(["--estimate-poly-a"])

    command.extend(["sup", reads])
    
    print(command)
    command_str = " ".join(command)  
    print(f">>> {command_str}\n")
    
    # launch job instantly. WILL NEED TO LAUNCH A JOB.
    subprocess.run(command, check=True)

    #saving(toml_config, "dorado")


if __name__ == "__main__":
    main()
