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


def main():
    parser = argparse.ArgumentParser(
        prog="LRS compare samples",
        description="Pipeline to compare LRS samples.",
    )

    parser.add_argument("config", type=str, help="Project config file.")

    args = parser.parse_args()

    # Loading initial TOML config
    with open(args.config, "r") as f:
        toml_config = toml.load(f)

    current_directory = os.getcwd()
    username = os.environ.get("USER")
    name = toml_config["general"]["comparison_name"]
    metadata = toml_config["general"]["metadata"]

    # Make output directory
    scripts = f"{current_directory}/scripts"
    if not os.path.exists(scripts):
        os.makedirs(scripts)

    scratch = f"/lustre10/scratch/{username}/{name}/"
    if not os.path.exists(scratch):
        os.makedirs(scratch)

    # Get functions to run
    function_queue = []

    if "degs" in toml_config["general"]["analyses"]:
        function_queue.append(flair_diffexp)

    if "as" in toml_config["general"]["analyses"]:
        df = pd.read_csv(metadata, sep="\t", header=0)
        make_manifest_combine(df)
        make_manifest_quantify(df, toml_config)

        function_queue.append(flair_diffsplice)

    if "apa" in toml_config["general"]["analyses"]:
        function_queue.append(isoquant_apalord)

    if "dmrs" in toml_config["general"]["analyses"]:
        function_queue.append(modkit)

    # Calling the functions
    for func in function_queue:
        func(toml_config)


def get_reference(ref):
    # Reference files location
    path = TOOL_PATH + "references/gencode/"
    reference: {}  # type: ignore
    match ref:
        case "grch38":
            reference = {
                "fasta": path + "GRCh38_p14/GRCh38.primary_assembly.genome.fa",
                "gtf": path + "GRCh38_p14/gencode.v48.primary_assembly.annotation.gtf",
                "chrom_size": path
                + "GRCh38_p14/GRCh38.primary_assembly.genome.chromSizes.txt",
                # Other genomes can be added later
            }
    return reference


def read_metadata(toml_config):
    metadata = toml_config["general"]["metadata"]
    df = pd.read_csv(metadata, sep="\t", header=0)

    return df


def make_manifest_combine(df):
    current_directory = os.getcwd()

    df["isoform"] = "isoform"
    df["bed"] = (
        df["project_path"]
        + "/results/flair/"
        + df["samples"]
        + "/"
        + df["samples"]
        + ".isoforms.bed"
    )
    df["fa"] = (
        df["project_path"]
        + "/results/flair/"
        + df["samples"]
        + "/"
        + df["samples"]
        + ".isoforms.fa"
    )
    df["readmap"] = (
        df["project_path"]
        + "/results/flair/"
        + df["samples"]
        + "/"
        + df["samples"]
        + ".isoform.read.map.txt"
    )

    df = df[["samples", "isoform", "bed", "fa", "readmap"]]
    df.to_csv(
        f"{current_directory}/manifest_combine.txt", header=False, sep="\t", index=False
    )


def make_manifest_quantify(df, toml_config):
    current_directory = os.getcwd()
    name = toml_config["general"]["comparison_name"]
    username = os.environ.get("USER")
    scratch = f"/lustre10/scratch/{username}/{name}/flair_diff"

    df["batch"] = (
        df["project_path"]
        .str.rstrip("/")
        .str.split("/")
        .str[-2]
        .str.split("_", n=1)
        .str[1]
    )
    df["fastq"] = scratch + "/fastq/" + df["samples"] + ".fastq"

    df = df[["samples", "phenotype", "batch", "fastq"]]
    df.to_csv(
        f"{current_directory}/manifest_quantify.txt",
        header=False,
        sep="\t",
        index=False,
    )


def flair_diffexp(toml_config):
    print()


def flair_diffsplice(toml_config):
    tool = "flair_diffsplice"
    current_directory = os.getcwd()
    name = toml_config["general"]["comparison_name"]
    email = toml_config["general"]["email"]
    genome = (
        TOOL_PATH + "references/gencode/GRCh38_p14/GRCh38.primary_assembly.genome.fa"
    )
    gtf = (
        TOOL_PATH
        + "references/gencode/GRCh38_p14/gencode.v48.primary_assembly.annotation.gtf"
    )
    conditionA = toml_config["general"]["conditionA"]
    conditionB = toml_config["general"]["conditionB"]

    # Get samples
    df = read_metadata(toml_config)
    samples = df["samples"].tolist()
    str_samples = " ".join(samples)

    job = current_directory + "/scripts/" + tool + ".slurm"
    with open(
        TOOL_PATH
        + "main_pipelines/long-read/LongReadSequencingONT/compare_samples/template_flair_diff.txt",
        "r",
    ) as f:
        slurm = f.read()
        slurm_filled = slurm.format(
            name,
            email,
            genome,
            gtf,
            current_directory,
            str_samples,
            toml_config["general"]["metadata"],
            conditionA,
            conditionB,
        )

        with open(job, "w") as o:
            o.write(slurm_filled)


def isoquant_apalord(toml_config):
    tool = "isoquant_apalord"
    current_directory = os.getcwd()
    name = toml_config["general"]["comparison_name"]
    email = toml_config["general"]["email"]
    genome = (
        TOOL_PATH + "references/gencode/GRCh38_p14/GRCh38.primary_assembly.genome.fa"
    )
    gtf = (
        TOOL_PATH
        + "references/gencode/GRCh38_p14/gencode.v48.primary_assembly.annotation.gtf"
    )
    conditionA = toml_config["general"]["conditionA"]
    conditionB = toml_config["general"]["conditionB"]

    # Get samples
    df = read_metadata(toml_config)

    samples = df["samples"].tolist()
    str_samples = " ".join(samples)

    cases = df[df["phenotype"] == conditionA]
    cases_list = cases["samples"].tolist()
    str_cases = ",".join(cases_list)

    ctrls = df[df["phenotype"] == conditionB]
    ctrls_list = ctrls["samples"].tolist()
    str_ctrls = ",".join(ctrls_list)

    # Run tool
    job = current_directory + "/scripts/" + tool + ".slurm"
    with open(
        TOOL_PATH
        + "main_pipelines/long-read/LongReadSequencingONT/compare_samples/template_isoquant_apalord.txt",
        "r",
    ) as f:
        slurm = f.read()
        slurm_filled = slurm.format(
            name,
            email,
            genome,
            gtf,
            str_samples,
            toml_config["general"]["metadata"],
            conditionA,
            conditionB,
            str_cases,
            str_ctrls,
        )

        with open(job, "w") as o:
            o.write(slurm_filled)


def modkit(toml_config):
    tool = "modkit_dmrs"
    current_directory = os.getcwd()
    name = toml_config["general"]["comparison_name"]
    email = toml_config["general"]["email"]
    genome = (
        TOOL_PATH + "references/gencode/GRCh38_p14/GRCh38.primary_assembly.genome.fa"
    )
    conditionA = toml_config["general"]["conditionA"]
    conditionB = toml_config["general"]["conditionB"]

    # Get samples
    df = read_metadata(toml_config)
    samples = df["samples"].tolist()
    str_samples = " ".join(samples)

    result_list = []
    case = df[df["phenotype"] == conditionA]
    ctrl = df[df["phenotype"] == conditionB]

    for row in ctrl.itertuples():
        line = (
            f"-a {row.project_path}/results/epi2me/{row.samples}/{row.samples}.wf_mods.1.bedmethyl.gz "
            f"-a {row.project_path}/results/epi2me/{row.samples}/{row.samples}.wf_mods.2.bedmethyl.gz"
        )
        result_list.append(line)

    for row in case.itertuples():
        line = (
            f"-b {row.project_path}/results/epi2me/{row.samples}/{row.samples}.wf_mods.1.bedmethyl.gz "
            f"-b {row.project_path}/results/epi2me/{row.samples}/{row.samples}.wf_mods.2.bedmethyl.gz"
        )
        result_list.append(line)

    print(result_list)
    pairs = " ".join(result_list)
    print(pairs)

    job = current_directory + "/scripts/" + tool + ".slurm"
    with open(
        TOOL_PATH
        + "main_pipelines/long-read/LongReadSequencingONT/compare_samples/template_modkit_dmrs.txt",
        "r",
    ) as f:
        slurm = f.read()
        slurm_filled = slurm.format(
            name,
            email,
            genome,
            toml_config["general"]["metadata"],
            str_samples,
            pairs,
        )

        with open(job, "w") as o:
            o.write(slurm_filled)


# Launches main function
if __name__ == "__main__":
    main()
