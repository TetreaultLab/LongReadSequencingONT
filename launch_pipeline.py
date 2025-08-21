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

    ## config and sample sheet
    if args.config != "config_final.toml":
        # Create final TOML config
        toml_config = create_config_final(args.config)

        # Create sample_sheet.csv
        create_sample_sheet(toml_config)

        print("\n\n\n!!! WARNING !!!\nIf you to change the parameters: Press CTRL+C now!\nModify config_final.toml and launch_pipeline with that config file.\n\nOtherwise it will run with default parameters.\n\n")
        time.sleep(10)

    else:
        toml_config = toml_config_initial
    
    # Get steps done
    f = open(output + "/scripts/steps_done.txt", "a")
    f.close()
    
    done = []
    with open(output + "/scripts/steps_done.txt", "r") as f:
        for line in f:
            done.append(line.strip())
    
    
    function_queue = []
    # Setting up list of steps
    # Base calling
    if "main_pipeline" not in done:
        function_queue.append(main_pipeline)

    #if "qc" not in done:
        #function_queue.append(qc)
    
    # SNP calling
    # if "SNP" in toml_config["general"]["analysis"]:
    #     if "clair3" not in done:
    #         function_queue.append(clair3)

    #     if "whatshap" not in done: # phasing
    #         function_queue.append(whatshap)

    # # SV calling
    # if "SV" in toml_config["general"]["analysis"]:
    #     if "sniffles2" not in done:
    #         function_queue.append(sniffles2)

    # Other tools ...


    # Create main.sh
    with open(output + "/scripts/main.sh", "w") as f:
        f.write("#!/bin/sh\n")

    # Calling each steps
    for func in function_queue:
        func(toml_config)

    # Call main.sh
    #subprocess.run(["bash", output + "/scripts/main.sh"])


def create_config_final(filename):
    # Read initial config file
    with open(filename, "r") as f:
        toml_config = toml.load(f)

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

    
    # get all flowcell directory names
    fc_dir_names = []
    for entry in os.listdir(new_path):
        full_path = os.path.join(new_path, entry)
        if os.path.isdir(full_path) and entry.startswith("2"):
            fc_dir_names.append(entry)
    toml_config['general']['fc_dir_names'] = fc_dir_names


    # Making directory structure in project
    directories = ["scripts", "alignments", "results", "qc"]
    for d in directories:
        if not os.path.exists(new_path + "/" + d):
            os.makedirs(new_path + "/" + d)

    # Making directory structure in flowcells subdirectories
    directories = ["main_reports", "reads", "alignments"]
    for flow in fc_dir_names:
        for d in directories:
            if not os.path.exists(new_path + "/" + flow + "/" + d):
                os.makedirs(new_path + "/" + flow + "/" + d)

        # Move main reports files to corresponding directory
        reports = ["barcode_alignment", "final_summary", "pore_activity", "report", "sample_sheet", "sequencing_summary", "throughput"]
        for r in reports:
            f = glob(os.path.join(new_path + "/" + flow + "/", r + "*"))
            if f != []: # if the file exist in the main directory
                for i in range(0, len(f)): # if one or more file starts with the report name
                    f2 = f[i]
                    name = f2.split("/")[-1]
                    os.rename(f2, os.path.join(new_path + "/" + flow + "/main_reports", name)) # move to main_reports

        # Move pod5 to reads directory
        if not os.path.exists(new_path + "/" + flow + "/reads/pod5"):
            os.rename(os.path.join(new_path, flow, "pod5"), os.path.join(new_path, flow, "reads", "pod5"))
        if os.path.exists(new_path + "/" + flow + "/pod5_skip"):
            for file in os.listdir(os.path.join(new_path, flow, "pod5_skip")):
                os.rename(os.path.join(new_path, flow, "pod5_skip", file), os.path.join(new_path, flow, "reads", "pod5", file))
            

    # Add Dorado options
    ## general options
    toml_config["dorado"] = {}
    toml_config['dorado']['min_q_score'] = 5
    toml_config['dorado']['barcode_both_ends'] = "False"
    toml_config['dorado']['model'] = "dna_r10.4.1_e8.2_400bps_sup@v5.2.0"

    ## parameters depending on sequencing type
    if seq_type == "WGS":
        toml_config['dorado']['mm2_opts'] = '"-x lr:hq"'
        if kit not in ["SQK-RBK114-24", "SQK-NBD114-24", "SQK-LSK114"]:
            raise Exception("Error: Wrong Kit for WGS. Options are SQK-RBK114-24, SQK-NBD114-24, SQK-LSK114")

    if seq_type == "RNA":
        toml_config['dorado']['mm2_opts'] = '"-x splice:hq"'
        if kit not in ["SQK-PCB114-24"]:
            raise Exception("Error: Wrong Kit for Whole Transcriptome. Options are SQK-PCB114-24")

    if methylation_status:
        toml_config['dorado']['modified_bases'] = "dna_r10.4.1_e8.2_400bps_sup@v5.2.0_5mC_5hmC@v1"
        toml_config['dorado']['modified_bases_threshold'] = 0.05

    if seq_type == "Targeted":
        toml_config['dorado']['mm2_opts'] = '"-x splice --junc-bed anno.bed12"'
        if kit not in ["SQK-NBD114-24"]:
            raise Exception("Error: Wrong Kit for Targeted Sequencing. Options are SQK-NBD114-24")


    # Add next tool options


    # Save new config file
    with open(toml_config["general"]["project_path"] + '/scripts/config_final.toml', 'w') as f:
        toml.dump(toml_config, f)

    return toml_config


def create_sample_sheet(toml_config):
    path = toml_config["general"]["project_path"]
    #rm_prefix = path.replace('/lustre09/project/6019267/shared/projects/Nanopore_Dock/', '')
    rm_prefix = path.replace('/lustre09/project/6019267/shared/tools/main_pipelines/long-read/', '')
    path_list = rm_prefix.split("/")
    project_name_date = path_list[0].split("_", 1)
    project_name = project_name_date[1]
    kit = toml_config["general"]["kit"]
    samples = toml_config["general"]["samples"]
    conditions = toml_config["general"]["conditions"]

    flowcells = toml_config["general"]["fc_dir_names"]
    barcode = toml_config["general"]["barcode"]

    for flowcell in flowcells:
        flow_cell_id_list = flowcell.split("_")
        flow_cell_id = flow_cell_id_list[3] + "_" + flow_cell_id_list[4]

        if type(barcode) is list:
            d = {'flow_cell_id': flow_cell_id, 
                'experiment_id': project_name, 
                'kit': kit, 
                'alias': samples, 
                'type': conditions, 
                'barcode': barcode}
        else :
            d = {'flow_cell_id': flow_cell_id, 
                'experiment_id': project_name, 
                'kit': kit, 
                'alias': samples, 
                'type': conditions, 
                'barcode': range(barcode, barcode + len(samples))}

        df = pd.DataFrame(data=d)
        df["barcode"] = "barcode" + df["barcode"].astype(int).astype(str).str.zfill(2)

        df.to_csv(path+"/scripts/" + flowcell + ".csv", sep=",", index=False)


def get_reference(ref):
    path = "/lustre09/project/6019267/shared/tools/references/gencode/"
    reference: {}  # type: ignore
    match ref:
        case "grch38":
            reference = {
                "fasta": path + "GRCh38_p14/GRCh38.primary_assembly.genome.fa",
                "gtf": path + "GRCh38_p14/gencode.v48.primary_assembly.annotation.gtf",
                "chrom_size": path + "GRCh38_p14/GRCh38.primary_assembly.genome.chromSizes.txt"
            }
    return reference


def create_script(tool, cores, memory, time, output, email, command, flowcell):
    if flowcell != "":
        job = output + "/scripts/" + tool + "_" + flowcell + ".slurm"

        with open("/lustre09/project/6019267/shared/tools/main_pipelines/long-read/LongReadSequencingONT/sbatch_template.txt", "r") as f:
            slurm = f.read()
            if tool == "dorado_basecaller":
                slurm_filled = slurm.format(cores, "#SBATCH --gpus=h100:1", memory, time, tool, flowcell, "log", "log", "def", email)
            
            else :
                slurm_filled = slurm.format(cores, "", memory, time, tool, flowcell, "log", "log", "rrg", email)
                slurm_filled += "module load StdEnv/2023 apptainer samtools\n"

            slurm_filled += "\n#\n### Calling " + tool + " - " + flowcell + "\n#\n"
            slurm_filled += command
            slurm_filled += "\n\n"
            # slurm_filled += 'echo "' + tool + '" >> ' + output + '/scripts/steps_done.txt'

            with open(job, "w") as o:
                o.write(slurm_filled)
            
                return job
    else :
        job = output + "/scripts/" + tool + ".slurm"

        with open("/lustre09/project/6019267/shared/tools/main_pipelines/long-read/LongReadSequencingONT/sbatch_template.txt", "r") as f:
            slurm = f.read()
            if tool == "dorado_demux":
                slurm_filled = slurm.format(cores, "", memory, time, tool, "run", "out", "log", "rrg", email)

            else :
                slurm_filled = slurm.format(cores, "", memory, time, tool, "run", "log", "log", "rrg", email)

            slurm_filled += "module load StdEnv/2023 samtools\n"
            slurm_filled += "source /lustre09/project/6019267/shared/tools/main_pipelines/long-read/launch_pipeline_env/bin/activate"

            slurm_filled += "\n#\n### Calling " + tool + "\n#\n"
            slurm_filled += command
            slurm_filled += "\n\n"
            # slurm_filled += 'echo "' + tool + '" >> ' + output + '/scripts/steps_done.txt'

            with open(job, "w") as o:
                o.write(slurm_filled)
            
                return job

def format_time(hours):
    days = int(hours // 24)
    remaining_hours = int(hours % 24)
    minutes = int((hours % 1) * 60)

    # Format as DD-HH:MM
    formatted_time = f"{days:02d}-{remaining_hours:02d}:{minutes:02d}"
    return(formatted_time)


def main_pipeline(toml_config):
    print("main_pipeline")
    output = toml_config["general"]["project_path"]
    email = toml_config["general"]["email"]
    genome = get_reference(toml_config["general"]["reference"])["fasta"]

    flowcells = toml_config["general"]["fc_dir_names"]
    codes = []
    for flowcell in flowcells:
        reads = output + "/" + flowcell + "/reads/pod5"
        final = output + "/" + flowcell + "/alignments/"
        bam_dorado = final + flowcell + ".bam"
        
        code = flowcell.split('_')[-1]
        codes.append(code)

        # BASECALLER
        tool = "dorado_basecaller"
        print(tool)
        cores = "8"
        memory = "64"

        # Get reads size 
        cmd = ["du", "-sh", "--apparent-size", "--block-size", "G", reads]
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        size_str = result.stdout.split()[0].rstrip('G')

        hours = int(size_str) * 0.02
        formatted_time = format_time(hours)

        command = ["/lustre09/project/6019267/shared/tools/main_pipelines/long-read/dorado-1.0.0-linux-x64/bin/dorado", "basecaller", "-v", "--device", "cuda:0", "--emit-moves", "--min-qscore", str(toml_config["dorado"]["min_q_score"]), "--reference", genome, "--sample-sheet", output + "/scripts/" + flowcell + ".csv", "--no-trim", "--kit-name", toml_config["general"]["kit"], "--mm2-opts", toml_config["dorado"]["mm2_opts"]]
        
        if toml_config["dorado"]["barcode_both_ends"] in ["true", "True", "yes", "Yes"]:
            command.extend(["--barcode-both-ends"])

        if "methylation" in toml_config["general"]["analysis"]:
            command.extend(["--modified-bases-models", "/lustre09/project/6019267/shared/tools/main_pipelines/long-read/dorado_models/" + toml_config["dorado"]["modified_bases"], "--modified-bases-threshold", str(toml_config["dorado"]["modified_bases_threshold"])])

        if "polya" in toml_config["general"]["analysis"]:
            command.extend(["--estimate-poly-a"])

        model = "/lustre09/project/6019267/shared/tools/main_pipelines/long-read/dorado_models/" + toml_config['dorado']['model']
        command.extend([model, reads, ">", bam_dorado])
        
        command_str = " ".join(command)

        job = create_script(tool, cores, memory, formatted_time, output, email, command_str, flowcell)

        var_name_bc = f"basecall_{code}"

        # Add slurm job to main.sh
        with open(output + "/scripts/main.sh", "a") as f:
            f.write("# Flowcell : " + flowcell + "\n")
            f.write(f'{var_name_bc}=$(sbatch --parsable ' + job + ")\n")

        # DEMUX
        tool2 = "dorado_demux"
        print(tool2)

        cores2 = "4"   
        memory2 = "24"

        hours2 = int(size_str) * 0.01
        formatted_time2 = format_time(hours2)
            
        command2 = ["/lustre09/project/6019267/shared/tools/main_pipelines/long-read/dorado-1.1.0-linux-x64/bin/dorado", "demux", "-vv", "--threads", cores2, "--no-trim", "--output-dir", final, "--no-classify", bam_dorado, "\n\n"]
        command_str2 = " ".join(command2)

        # Create slurm job
        job2 = create_script(tool2, cores2, memory2, formatted_time2, output, email, command_str2, flowcell)
        
        var_name = f"demux_{code}"
        
        # Add slurm job to main.sh
        with open(output + "/scripts/main.sh", "a") as f:
            f.write(f'{var_name}=$(sbatch --parsable --dependency=afterok:${var_name_bc} ' + job2 + ')\n\n')

    # Samtools
    print("samtools")
    dependencies = ":".join([f"$demux_{code}" for code in codes])
    command3 = ["python", "-u", "/lustre09/project/6019267/shared/tools/main_pipelines/long-read/LongReadSequencingONT/rename_bam.py", toml_config["general"]["project_path"] + '/scripts/config_final.toml']
    command_str3 = " ".join(command3)
    time3 = "00-11:00"
    job3 = create_script("samtools", "8", "32", time3, output, email, command_str3, "")

    with open(output + "/scripts/main.sh", "a") as f:
        f.write("# Rename, merge, sort and index bams")
        f.write(f"\nsamtools=$(sbatch --parsable --dependency=afterok:{dependencies} {job3})\n")


    # QC
    print("QC")
    command_str4 = ""
    for name in toml_config["general"]["samples"]:
        command4 = ["apptainer", "run", "/lustre09/project/6019267/shared/tools/main_pipelines/long-read/image_longreadsum.sif", "bam", "--threads", "8", "--log", output + "/qc/longreadsum_" + name + ".log", "--ref", genome, "-Q", '"' + name + '_"', "-i", output + "/alignments/" + name + "_sorted.bam", "-o", output + "/qc"]

        if "methylation" in toml_config["general"]["analysis"]:
            command4.extend(["--mod"])
        
        command4.extend(["\n\n"])
        command_str4 += " ".join(command4)

    time4 = "00-11:00"

    job4 = create_script("LongReadSum", "8", "64", time4, output, email, command_str4, "")

    with open(output + "/scripts/main.sh", "a") as f:
        f.write("# QC\n")
        f.write(f'\nsbatch --dependency=afterok:$samtools {job4})\n')



def qc(toml_config):
    tool = "qc"
    output = toml_config["general"]["project_path"]
    df = pd.read_csv(output + "/scripts/samples.csv", header=0)
    df["name"] = df["alias"] + "_" + df["barcode"]
    fasta = get_reference(toml_config["general"]["reference"])["fasta"]
    
    threads = "8"
    memory = "200"
    time = "00-11:59"
    email = toml_config["general"]["email"]

    command_pod5 = ["apptainer", "run", "/lustre03/project/6019267/shared/tools/PIPELINES/LongReadSequencing/image_longreadsum.sif", "pod5", "--threads", threads, "--log", output + "/qc/longreadsum_pod5.log", "-Q", '"' + flowcell + '"', "-P", '"' + output + "/reads/pod5/*.pod5\"", "-o", output + "/qc", "--basecalls", output + "/alignments/" + flowcell + ".bam", "\n\n"]
    command_str1 = " ".join(command_pod5)
    
    command_str2 = ""
    for name in df["name"]:
        command = ["apptainer", "run", "/lustre03/project/6019267/shared/tools/PIPELINES/LongReadSequencing/image_longreadsum.sif", "bam", "--threads", threads, "--log", output + "/qc/longreadsum_" + name + ".log", "--ref", fasta, "-Q", '"' + name + '_"', "-i", output + "/alignments/" + name + ".bam", "-o", output + "/qc"]

        if "methylation" in toml_config["general"]["analysis"]:
            command.extend(["--mod"])
        
        command.extend(["\n\n"])
        command_str2 += " ".join(command)

    command_str = command_str1 + command_str2
    
    # Create slurm job
    job = create_script(tool, threads, memory, time, output, email, command_str)

    done = []
    with open(output + "/scripts/steps_done.txt", "r") as f:
        for line in f:
            done.append(line.strip())
    
    if "dorado_demux" not in done:
        with open(output + "/scripts/main.sh", "a") as f:
            f.write("sbatch --dependency=afterok:$dorado " + job + "\n")
    else:
        with open(output + "/scripts/main.sh", "a") as f:
            f.write("sbatch " + job + "\n")


def clair3(toml_config):
    tool = "clair3"

    output = toml_config["general"]["project_path"]
    ref = get_reference(toml_config["general"]["reference"])["fasta"] # same ref for RNA + DNA ?
    bam = output + "/alignments/" + flowcell + "_sorted.bam"
    model = "/lustre03/project/6019267/shared/tools/PIPELINES/LongReadSequencing/Clair3/models/r1041_e82_400bps_sup_v420" # https://www.bio8.cs.hku.hk/clair3/clair3_models/
    platform_rna = "ont_dorado_drna004" # Possible options: {ont_dorado_drna004, ont_guppy_drna002, ont_guppy_cdna, hifi_sequel2_pbmm2, hifi_sequel2_minimap2, hifi_mas_pbmm2, hifi_sequel2_minimap2}.
    platform_dna = "ont"

    threads = "8"
    memory = "32"
    time = "00-23:59"
    email = toml_config["general"]["email"]

    if toml_config["general"]["seq_type"] == "RNA":
        command = ["apptainer", "run", "-C", "-W", "${SLURM_TMPDIR}", "/lustre03/project/6019267/shared/tools/PIPELINES/LongReadSequencing/image_" + tool + ".sif ", "/lustre03/project/6019267/shared/tools/PIPELINES/LongReadSequencing/Clair3-RNA/run_clair3_rna", "--bam_fn", bam, "--ref_fn", ref, "--threads", threads, "--platform", platform_rna, "--output_dir", output + "/results/", "\n\n"]
        # add --enable_phasing_model ?
    else:
        command = ["apptainer", "run", "-C", "-W", "${SLURM_TMPDIR}", "/lustre03/project/6019267/shared/tools/PIPELINES/LongReadSequencing/image_" + tool + ".sif ", "run_clair3.sh", "-b", bam, "-f", ref, "-m", model, "-t", threads, "-p", platform_dna, "-o", output + "/results/", "\n\n"]

    command2 = ["mv", output + "/results/output.vcf.gz", output + "/results/" + flowcell + ".vcf.gz", "\n\n"]
    command3 = ["mv", output + "/results/output.vcf.gz.tbi", output + "/results/" + flowcell + ".vcf.gz.tbi"]
    
    command_str1 = " ".join(command)
    command_str2 = " ".join(command2)
    command_str3 = " ".join(command3)

    command_str = command_str1 + command_str2 + command_str3

    # Create slurm job
    job = create_script(tool, threads, memory, time, output, email, command_str)
    

    # Add slurm job to main.sh
    done = []
    with open(output + "/scripts/steps_done.txt", "r") as f:
        for line in f:
            done.append(line.strip())

    if "dorado" not in done:
        with open(output + "/scripts/main.sh", "a") as f:
            f.write("sbatch --dependency=afterok:$dorado " + job + "\n")
    else:
        with open(output + "/scripts/main.sh", "a") as f:
            f.write("sbatch " + job + "\n")


def whatshap(toml_config):
    tool = "whatshap"

    output = toml_config["general"]["project_path"]
    output_vcf = output + "/results/" + flowcell + "phased.vcf"
    input_vcf = output + "/results/" + flowcell + ".vcf.gz"
    bam = output + "/alignments/" + flowcell + "_sorted.bam"
    
    ref = get_reference(toml_config["general"]["reference"])["fasta"]

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


def sniffles2(toml_config):
    tool = "sniffles2"

    output = toml_config["general"]["project_path"]
    bam = "/home/shared/data/2024-10-16_Lapiana_n17/no_sample_id/20241016_1653_X2_FAV26227_d404da0e/alignment/minimap2_sup/B1540_sorted.bam"
    vcf = output + "/sniffles2_SV.vcf"
    ref = get_reference(toml_config["general"]["reference"])["fasta"]

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


if __name__ == "__main__":
    main()
