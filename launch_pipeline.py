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
    
    # SNP calling
    if "clair3" not in done and "clair3_rna" not in done:
        print(">>> Variant calling - SNP: Clair3 (?)")
        function_queue.append(clair3)

    # Phasing
    if "whatshap" not in done:
        print(">>> Variant phasing: WhatsHap (?)")
        function_queue.append(whatshap)

    # SV calling
    # if "sniffles2" not in done:
    #     print(">>> Variant calling - SV: Sniffles2 (?)")
    #     function_queue.append(sniffles2)

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


def get_reference(ref, tool):
    path = "/lustre03/project/6019267/shared/tools/database_files/hg38/"
    reference: {}  # type: ignore
    match ref:
        case "grch38":
            reference = {
                "fasta": path + "gencode.v38.p13.genome.fa",
                "gtf": path + "gencode.v38.annotation.gtf",
            }

    return reference


def create_script(tool, cores, memory, time, output, email, command):
    project_name = "" # output.split() # TO-DO
    job = output + "/scripts/" + tool + ".slurm"

    with open("/lustre03/project/6019267/shared/tools/PIPELINES/LongReadSequencing/LongReadSequencingONT/sbatch_template.txt", "r") as f:
        slurm = f.read()
        slurm_filled = slurm.format(cores, memory, time, tool, project_name, email)
        
        slurm_filled += "module load StdEnv/2023 dorado/0.8.3 apptainer"

        slurm_filled += "\n#\n### Calling " + tool + "\n#\n"
        slurm_filled += "apptainer run -C -W ${SLURM_TMPDIR} --nv -B /project -B /scratch " + command
        slurm_filled += "\n"

        with open(job, "w") as o:
            o.write(slurm_filled)
        
        return job


def dorado(toml_config):
    tool = "dorado"
    title(tool)
    
    output = toml_config["general"]["project_path"]
    email = toml_config["general"]["email"]
    genome = get_reference(toml_config["general"]["reference"], tool)["fasta"]
	reads = "/lustre03/project/6019267/shared/projects/Nanopore_Dock/pod5"
    #reads = "/home/shared/data/2025-01-15_FXN-Batch4/FRDA14_21-UTMAB-06_2/20250115_2147_P2S-02441-B_PBA20836_7fce705b/pod5/PBA20836_7fce705b_9c89ba7b_66.pod5"
    
    cores = 8
    memory = 32
    time = "02-23:59"

    command = ["dorado", "basecaller", "--verbose", "--device", "cuda:0", "--min-qscore", str(toml_config["dorado"]["min_q_score"]), "-o", output, "--reference", genome, "--sample-sheet", output + "/" + toml_config["dorado"]["sample_sheet"], "--trim", toml_config["dorado"]["trim"], "--kit-name", toml_config["general"]["kit"], "--mm2-opts", toml_config["dorado"]["mm2_opts"]]
    
    if toml_config["dorado"]["barcode_both_ends"] in ["true", "True", "yes", "Yes"]:
        command.extend(["--barcode-both-ends"])

    if "methylation" in toml_config["general"]["analysis"]:
        command.extend(["--modified-bases", toml_config["dorado"]["modified_bases"], "--modified-bases-threshold", str(toml_config["dorado"]["modified_bases_threshold"])])

    if "polya" in toml_config["general"]["analysis"]:
        command.extend(["--estimate-poly-a"])

    command.extend(["sup", reads])
    
    command_str = " ".join(command)  
    print(f">>> {command_str}\n")
    
    # Create slurm job
    job = create_script(tool, cores, memory, time, output, email, command_str)
    
    # Launch slurm job
    subprocess.run(["bash", job], check=True) # put sbatch instead of bash when on beluga
    
    # Mark tool as done
    saving(toml_config, tool)


def clair3(toml_config):
    if toml_config["general"]["seq_type"] == "RNA":
        tool = "clair3_rna"
    else:
        tool = "clair3"
    
    title(tool)

    output = toml_config["general"]["project_path"]
    ref = get_reference(toml_config["general"]["reference"], tool)["fasta"] # same ref for RNA + DNA ?
    bam = output + ".bam" # complete with file name from dorado
    #bam = "/home/shared/data/2024-10-16_Lapiana_n17/no_sample_id/20241016_1653_X2_FAV26227_d404da0e/alignment/minimap2_sup/B1540_sorted.bam"
    model = "/home/shared/tools/clair3/Clair3/models/r1041_e82_400bps_sup_v420" # https://www.bio8.cs.hku.hk/clair3/clair3_models/
    platform_rna = "ont_dorado_drna004" # Possible options: {ont_dorado_drna004, ont_guppy_drna002, ont_guppy_cdna, hifi_sequel2_pbmm2, hifi_sequel2_minimap2, hifi_mas_pbmm2, hifi_sequel2_minimap2}.
    platform_dna = "ont"

    threads = "8"
    memory = "32"
    time = "00-23:59"
    email = toml_config["general"]["email"]

    if tool == "clair3_rna":
        command = ["run_clair3_rna", "--bam_fn", bam, "--ref_fn", ref, "--threads", threads, "--platform", platform_rna, "--output_dir", output]
        # add --enable_phasing_model ?
    else:
        command = ["run_clair3.sh", "-b", bam, "-f", ref, "-m", model, "-t", threads, "-p", platform_dna, "-o", output]

    command_str = " ".join(command)  
    print(f">>> {command_str}\n")

    # Create slurm job
    job = create_script(tool, threads, memory, time, output, email, command_str)
    
    # Launch slurm job
    subprocess.run(["bash", job], check=True) # put sbatch instead of bash when on beluga
    
    # Mark tool as done
    saving(toml_config, tool)


def whatshap(toml_config):
    tool = "whatshap"
    title(tool)

    output = toml_config["general"]["project_path"] 
    output_vcf = output + "/phased.vcf"
    input_vcf = output + "/merge_output.vcf.gz"
	bam = output + "_sorted.bam" # complete with file name from dorado
    bam = "/home/shared/data/2024-10-16_Lapiana_n17/no_sample_id/20241016_1653_X2_FAV26227_d404da0e/alignment/minimap2_sup/B1540_sorted.bam"
    ref = get_reference(toml_config["general"]["reference"], tool)["fasta"]

    threads = "8"
    memory = "32"
    time = "00-23:59"
    email = toml_config["general"]["email"]

    command = ["whatshap", "phase", "--ignore-read-groups", "-o", output_vcf, "--reference", ref, input_vcf, bam]

    command_str = " ".join(command)  
    print(f">>> {command_str}\n")

    # Create slurm job
    job = create_script(tool, threads, memory, time, output, email, command_str)
    
    # Launch slurm job
    subprocess.run(["bash", job], check=True) # put sbatch instead of bash when on beluga
    
    # Mark tool as done
    saving(toml_config, tool)


def sniffles2(toml_config):
    tool = "sniffles2"
    title(tool)

    output = toml_config["general"]["project_path"]
    threads = "8"
    memory = "32"
    time = "00-23:59"
    email = toml_config["general"]["email"]

    # to-do
    command = []

    command_str = " ".join(command)  
    print(f">>> {command_str}\n")

    # Create slurm job
    job = create_script(tool, threads, memory, time, output, email, command_str)
    
    # Launch slurm job
    subprocess.run(["bash", job], check=True) # put sbatch instead of bash when on beluga
    
    # Mark tool as done
    #saving(toml_config, tool)

if __name__ == "__main__":
    main()
