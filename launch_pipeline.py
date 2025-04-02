import argparse
import toml
import time
from datetime import datetime
import subprocess
import sys
import pandas as pd
import numpy as np
import traceback
import os
from glob import glob

def main():
    start = get_time()

    # Start of pipeline
    start_str = ">>> LRS pipeline starting at {}.".format(
        start
    )
    print(
        "=" * len(start_str) + "\n" + start_str + "\n" + "=" * len(start_str),
        file=sys.stdout,
    )
    
    parser = argparse.ArgumentParser(
        prog="PipelineLong",
        description="LongReadSequencing pipeline.",
    )

    parser.add_argument(
        "config", type=str, help="Project config file, including path."
    )

    args = parser.parse_args()
    
    # Loading initial TOML config
    with open(args.config, "r") as f:
        toml_config_initial = toml.load(f)

    output = toml_config_initial["general"]["project_path"]
    
    # Create list of steps already done from the file steps_done.txt
    steps = open(output + "/steps_done.txt", "a").close()
    done = []
    with open(output + "/steps_done.txt", "r") as f:
        for line in f:
            done.append(line.strip())

    # Creat config and sample sheet
    if "Loading ENV" not in done:
        # Create final TOML config
        toml_config = create_config_final(args.config)

        # Create sample_sheet.csv
        create_sample_sheet(toml_config)

        # Open file for steps done
        steps = open(output + "/steps_done.txt", "a")
        steps.write("Loading ENV\n")
        steps.close()

        print("\n\n\n!!! WARNING !!!\nIf you to change the parameters: Press CTRL+C now!\nModify config_final.toml and launch_pipeline with that config file.\n\nOtherwise it will run with default parameters.\n\n")
        time.sleep(60)

    else:
        toml_config = toml_config_initial
    

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
    # if "clair3" not in done and "clair3_rna" not in done:
    #    print(">>> Variant calling - SNP: Clair3 (?)")
    #    function_queue.append(clair3)

    # Phasing
    # if "whatshap" not in done:
    #     print(">>> Variant phasing: WhatsHap (?)")
    #     function_queue.append(whatshap)

    # SV calling
    # if "sniffles2" not in done:
    #    print(">>> Variant calling - SV: Sniffles2 (?)")
    #    function_queue.append(sniffles2)

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

def create_config_final(filename):
    # Read initial config file
    with open(filename, "r") as f:
        toml_config = toml.load(f)

    ## Which analysis are possible for each seq-type. Make a if to fill the analysis field and check if some error are made in analysis was not empty
    # if analysis list is empty, add all available analysis
    if toml_config["general"]["analysis"] == []:
        toml_config["general"]["analysis"] = ["methylation", "splicing", "polya", "snp", "sv", "repeats"]

    # set kit to variable
    kit = toml_config['general']['kit']

    # check seq_type
    seq_type = toml_config["general"]["seq_type"]
    if (toml_config["general"]["seq_type"] == "Meth") or ("methylation" in toml_config["general"]["analysis"]):
        methylation_status = True
    else:
        methylation_status = False

    # if path has a "/" at the end, remove it 
    path = toml_config["general"]["project_path"]
    if path.endswith('/'):
        new_path = path[0:-1]
        toml_config["general"]["project_path"] = new_path
    else:
        new_path = path

    # Making directory structure
    directories = ["main_reports", "reads", "scripts", "alignments", "results", "qc"]
    for d in directories:
        if not os.path.exists(new_path + "/" + d):
            os.makedirs(new_path + "/" + d)

    # Move main reports files to corresponding directory
    reports = ["barcode_alignment", "final_summary", "pore_activity", "report", "sample_sheet", "sequencing_summary", "throughput"]
    for r in reports:
        f = glob(os.path.join(new_path, r + "*"))
        if f != []: # if the file exist in the main directory
            for i in range(0, len(f)): # if one or more file starts with the report name
                f2 = f[i]
                name = f2.split("/")[-1]
                os.rename(f2, os.path.join(new_path, "main_reports", name)) # move to main_reports

    # Add Dorado options
    ## general options
    toml_config["dorado"] = {}
    toml_config['dorado']['min_q_score'] = 10
    toml_config['dorado']['sample_sheet'] = "samples.csv"
    toml_config['dorado']['barcode_both_ends'] = True
    toml_config['dorado']['trim'] = "all"
    toml_config['dorado']['model'] = "dna_r10.4.1_e8.2_400bps_sup@v5.0.0"

    ## parameters depending on sequencing type
    if seq_type == "WGS":
        toml_config['dorado']['mm2_opts'] = "'-ax lr:hq'"
        if kit not in ["SQK-RBK114-24", "SQK-NBD114.24", "SQK-LSK114"]:
            raise Exception("Error: Wrong Kit for WGS. Options are SQK-RBK114-24, SQK-NBD114.24, SQK-LSK114")

    if seq_type == "RNA":
        toml_config['dorado']['mm2_opts'] = "'-ax splice:hq -uf'"
        if kit not in ["SQK-PCB114-24"]:
            raise Exception("Error: Wrong Kit for Whole Transcriptome. Options are SQK-PCB114-24")

    if methylation_status:
        toml_config['dorado']['modified_bases'] = "5mCG_5hmCG"
        toml_config['dorado']['modified_bases_threshold'] = 0.05

    if seq_type == "Targeted":
        toml_config['dorado']['mm2_opts'] = "'-ax splice --junc-bed anno.bed12'"
        if kit not in ["SQK-NBD114-24"]:
            raise Exception("Error: Wrong Kit for Targeted Sequencing. Options are SQK-NBD114-24")


    # Add next tool options


    # Save new config file
    with open(toml_config["general"]["project_path"] + '/config_final.toml', 'w') as f:
        toml.dump(toml_config, f)

    return toml_config

def create_sample_sheet(toml_config):
    path = toml_config["general"]["project_path"]
    rm_prefix = path.replace('/lustre03/project/6019267/shared/projects/Nanopore_Dock/', '')
    path_list = rm_prefix.split("/")
    project_name_date = path_list[0].split("_", 1)
    project_name = project_name_date[1]
    flow_cell = path_list[2]
    flow_cell_id_list = flow_cell.split("_")
    flow_cell_id = flow_cell_id_list[3] + "_" + flow_cell_id_list[4]
    kit = toml_config["general"]["kit"]
    samples = toml_config["general"]["samples"]
    conditions = toml_config["general"]["conditions"]
    barcode_initial = toml_config["general"]["barcode"]

    d = {'flow_cell_id': flow_cell_id, 
         'experiment_id': project_name, 
         'kit': kit, 
         'alias': samples, 
         'type': conditions, 
         'barcode': range(barcode_initial, barcode_initial + len(samples))}

    df = pd.DataFrame(data=d)

    df.to_csv(path+"/samples.csv", sep=",", index=False)



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
    job = output + "/scripts/" + tool + ".slurm"
    rm_prefix = output.replace('/lustre03/project/6019267/shared/projects/Nanopore_Dock/', '')
    path_list = rm_prefix.split("/")
    project_name = path_list[0]

    with open("/lustre03/project/6019267/shared/tools/PIPELINES/LongReadSequencing/LongReadSequencingONT/sbatch_template.txt", "r") as f:
        slurm = f.read()
        # TO-DO: def-tetreaum for dorado and rrg-tetreaum for the rest
        slurm_filled = slurm.format(cores, memory, time, tool, project_name, email)
        
        slurm_filled += "module load StdEnv/2023 dorado/0.8.3 apptainer"

        slurm_filled += "\n#\n### Calling " + tool + "\n#\n"
        if tool in ["clair3", "clair3_rna", "whatshap"]:
       	    slurm_filled += "apptainer run -C -W ${SLURM_TMPDIR} --nv -B /project -B /scratch /lustre03/project/6019267/shared/tools/PIPELINES/LongReadSequencing/image_" + tool + ".sif " + command
        else:
	        slurm_filled += command
        
        slurm_filled += "\n"

        with open(job, "w") as o:
            o.write(slurm_filled)
        
            return job


def dorado(toml_config):
    tool = "dorado"
    title(tool)
    
    output = toml_config["general"]["project_path"]
    reads = output + "/pod5"
    final = output + "/alignments/"
    email = toml_config["general"]["email"]
    genome = get_reference(toml_config["general"]["reference"], tool)["fasta"]
    
    cores = 1
    memory = 12
    time = "00-10:59"

    command = ["dorado", "basecaller", "--verbose", "--device", "cuda:auto", "--min-qscore", str(toml_config["dorado"]["min_q_score"]), "--output-dir", final, "--reference", genome, "--sample-sheet", output + "/" + toml_config["dorado"]["sample_sheet"], "--trim", toml_config["dorado"]["trim"], "--kit-name", toml_config["general"]["kit"], "--mm2-opts", toml_config["dorado"]["mm2_opts"]]
    
    if toml_config["dorado"]["barcode_both_ends"] in ["true", "True", "yes", "Yes"]:
        command.extend(["--barcode-both-ends"])

    if "methylation" in toml_config["general"]["analysis"]:
        command.extend(["--modified-bases", toml_config["dorado"]["modified_bases"], "--modified-bases-threshold", str(toml_config["dorado"]["modified_bases_threshold"])])

    if "polya" in toml_config["general"]["analysis"]:
        command.extend(["--estimate-poly-a"])

    command.extend(["/lustre03/project/6019267/shared/tools/PIPELINES/LongReadSequencing/dorado_models/" + toml_config['dorado']['model'], reads])
    
    command_str = " ".join(command)  
    print(f">>> {command_str}\n")
    
    # Create slurm job
    job = create_script(tool, cores, memory, time, output, email, command_str)
    
    # Launch slurm job
    subprocess.run(["sbatch", job], check=True) # put sbatch instead of bash when on beluga
    
    # Mark tool as done
    #saving(toml_config, tool)


def clair3(toml_config):
    if toml_config["general"]["seq_type"] == "RNA":
        tool = "clair3_rna"
    else:
        tool = "clair3"
    
    title(tool)

    output = toml_config["general"]["project_path"]
    ref = get_reference(toml_config["general"]["reference"], tool)["fasta"] # same ref for RNA + DNA ?
    bam = output + "/alignments/" # complete with file name from dorado
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
    # bam = "/home/shared/data/2024-10-16_Lapiana_n17/no_sample_id/20241016_1653_X2_FAV26227_d404da0e/alignment/minimap2_sup/B1540_sorted.bam"
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
    bam = "/home/shared/data/2024-10-16_Lapiana_n17/no_sample_id/20241016_1653_X2_FAV26227_d404da0e/alignment/minimap2_sup/B1540_sorted.bam"
    vcf = output + "/sniffles2_SV.vcf"
    ref = get_reference(toml_config["general"]["reference"], tool)["fasta"]

    threads = "8"
    memory = "32"
    time = "00-23:59"
    email = toml_config["general"]["email"]

    command = ["sniffles", "--threads", threads, "--reference", ref, "--input", bam, "--vcf", vcf]

    command_str = " ".join(command)  
    print(f">>> {command_str}\n")

    # Create slurm job
    job = create_script(tool, threads, memory, time, output, email, command_str)
    
    # Launch slurm job
    subprocess.run(["bash", job], check=True) # put sbatch instead of bash when on beluga
    
    # Mark tool as done
    saving(toml_config, tool)


if __name__ == "__main__":
    main()
