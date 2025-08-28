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

# Global variables
TOOL_PATH = "/lustre09/project/6019267/shared/tools/"

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

    ## Create config and sample sheet
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
    steps_done_file = os.path.join(output, "scripts", "steps_done.txt")
    if os.path.exists(steps_done_file):
        with open(steps_done_file, "r") as f:
            done = [line.strip() for line in f if line.strip()]  # ignore empty lines

    function_queue = []
    # Setting up list of steps
    # Dorado Basecalling
    if "dorado_basecaller" not in done:
        function_queue.append(dorado_basecaller)

    # Dorado Basecalling
    if "dorado_demux" not in done:
        function_queue.append(dorado_demux)

    # Renaming bams and running samtools merge, sort and index
    if "samtools" not in done:
        function_queue.append(samtools)
    
    # Quality Control - LongReadSum
    if "longreadsum" not in done:
        function_queue.append(longreadsum)

    # Cleanup
    if "cleanup" not in done:
        function_queue.append(cleanup)


    # Create main.sh
    with open(output + "/scripts/main.sh", "w") as f:
        f.write("#!/bin/sh\n")

    # Calling each steps
    for func in function_queue:
        func(toml_config)

    # Call main.sh (Launch the pipeline)
    subprocess.run(["bash", output + "/scripts/main.sh"])


def create_config_final(filename):
    
    # Read initial config file
    with open(filename, "r") as f:
        toml_config = toml.load(f)

    # Set kit to variable
    kit = toml_config['general']['kit']

    # Check seq_type
    seq_type = toml_config["general"]["seq_type"]
    if (toml_config["general"]["seq_type"] == "Meth") or ("methylation" in toml_config["general"]["analysis"]):
        methylation_status = True
    else:
        methylation_status = False

    # Make sure path variable is correctly formatted (remove "/" at the end)
    path = toml_config["general"]["project_path"].rstrip("/")
    toml_config["general"]["project_path"] = path

    # Get all flowcell directory names
    fc_dir_names = []
    for entry in os.listdir(path):
        full_path = os.path.join(path, entry)
        # This assumes MinKNOW will not change flowcell naming convention
        if os.path.isdir(full_path) and entry.startswith("2"):
            fc_dir_names.append(entry)
    toml_config['general']['fc_dir_names'] = fc_dir_names


    # Making directory structure in project (Within each sample library)
    directories = ["scripts/logs", "alignments", "results", "qc"]
    for d in directories:
        if not os.path.exists(path + "/" + d):
            os.makedirs(path + "/" + d)

    # Making directory structure in flowcells subdirectories
    directories = ["main_reports", "reads", "alignments"]
    for flow in fc_dir_names:
        for d in directories:
            if not os.path.exists(path + "/" + flow + "/" + d):
                os.makedirs(path + "/" + flow + "/" + d)

        # Move main reports files to corresponding directory
        reports = ["barcode_alignment", "final_summary", "pore_activity", "report", "sample_sheet", "sequencing_summary", "throughput"]
        for r in reports:
            f = glob(os.path.join(path + "/" + flow + "/", r + "*"))
            if f != []: # if the file exist in the main directory
                for i in range(0, len(f)): # if one or more file starts with the report name
                    f2 = f[i]
                    name = f2.split("/")[-1]
                    os.rename(f2, os.path.join(path + "/" + flow + "/main_reports", name)) # move to main_reports

        # Move pod5 to reads directory
        if not os.path.exists(path + "/" + flow + "/reads/pod5"):
            os.rename(os.path.join(path, flow, "pod5"), os.path.join(path, flow, "reads", "pod5"))
        if os.path.exists(path + "/" + flow + "/pod5_skip"):
            for file in os.listdir(os.path.join(path, flow, "pod5_skip")):
                os.rename(os.path.join(path, flow, "pod5_skip", file), os.path.join(path, flow, "reads", "pod5", file))


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

    path = toml_config["general"]["project_path"].rstrip("/")
    parts = path.split(os.sep)

    if len(parts) < 2:
        raise ValueError(f"Project path '{path}' is too short to parse project/date info")

    # Take the *second-to-last* directory as project folder 
    project_folder = parts[-2]

    # Split date and project name
    try:
        project_date, project_name = project_folder.split("_", 1)
    except ValueError:
        raise ValueError(
            f"Project folder name '{project_folder}' does not contain an underscore "
            "to split date/project name"
        )

    # Extract necessary information from the user-defined configuration
    kit = toml_config["general"]["kit"]
    samples = toml_config["general"]["samples"]
    conditions = toml_config["general"]["conditions"]
    flowcells = toml_config["general"]["fc_dir_names"]
    barcode = toml_config["general"]["barcode"]

    # Iterate through each flowcell
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

        # Create the sample sheet and save to .csv
        df = pd.DataFrame(data=d)
        df["barcode"] = "barcode" + df["barcode"].astype(int).astype(str).str.zfill(2)

        df.to_csv(path+"/scripts/" + flowcell + ".csv", sep=",", index=False)


def get_reference(ref):

    # Reference files location
    path = TOOL_PATH + "references/gencode/"
    reference: {}  # type: ignore
    match ref:
        case "grch38":
            reference = {
                "fasta": path + "GRCh38_p14/GRCh38.primary_assembly.genome.fa",
                "gtf": path + "GRCh38_p14/gencode.v48.primary_assembly.annotation.gtf",
                "chrom_size": path + "GRCh38_p14/GRCh38.primary_assembly.genome.chromSizes.txt"
        # Other genomes can be added later
            }
    return reference


def create_script(tool, cores, memory, time, output, email, command, flowcell):

    # Enables creating a script per flowcell, or a single script if "" is added as the argument
    if flowcell != "":
        job = output + "/scripts/" + tool + "_" + flowcell + ".slurm"

        # Uses a slurm template for each job script
        with open(TOOL_PATH + "main_pipelines/long-read/LongReadSequencingONT/sbatch_template.txt", "r") as f:
            slurm = f.read()
            # Because dorado utilizes a GPU, it needs different options
            if tool == "dorado_basecaller":
                slurm_filled = slurm.format(cores, "#SBATCH --gpus=h100:1", memory, time, tool, flowcell, "log", "log", "def", email)
            
            else if tool == "dorado_demux":
                slurm_filled = slurm.format(cores, "", memory, time, tool, "run", "out", "log", "rrg", email)
            
            # Most tools will be CPU-dependent 
            else :
                slurm_filled = slurm.format(cores, "", memory, time, tool, flowcell, "log", "log", "rrg", email)
                slurm_filled += "module load StdEnv/2023 apptainer samtools\n"

            slurm_filled += "\n#\n### Calling " + tool + " - " + flowcell + "\n#\n"
            slurm_filled += command
            slurm_filled += "\n\n"
            # Pending completion to keep track of completed steps 
            #slurm_filled += 'echo "' + tool + '" >> ' + output + '/scripts/steps_done.txt'

            with open(job, "w") as o:
                o.write(slurm_filled)

                return job
            
    # This is for tools running on all the data at once
    else :
        job = output + "/scripts/" + tool + ".slurm"

        with open(TOOL_PATH + "main_pipelines/long-read/LongReadSequencingONT/sbatch_template.txt", "r") as f:
            slurm = f.read()
            # Because dorado has a different log format? (Not sure if obsolete as demux is flowcell-specific)
            # For all other tools 
            slurm_filled = slurm.format(cores, "", memory, time, tool, "run", "log", "log", "rrg", email)

            # Add enviroment loading commands
            slurm_filled += "module load StdEnv/2023 apptainer samtools\n"
            slurm_filled += "source " + TOOL_PATH + "main_pipelines/long-read/launch_pipeline_env/bin/activate"

            slurm_filled += "\n#\n### Calling " + tool + "\n#\n"
            slurm_filled += command
            slurm_filled += "\n\n"
            # Pending completion to keep track of completed steps
            # slurm_filled += 'echo "' + tool + '" >> ' + output + '/scripts/steps_done.txt'

            with open(job, "w") as o:
                o.write(slurm_filled)

                return job

def format_time(hours):

    # Job time allocation formatting
    if hours < 1:
        hours = 1.5
    days = int(hours // 24)
    remaining_hours = int(hours % 24)
    minutes = int((hours % 1) * 60)

    # Format as DD-HH:MM
    formatted_time = f"{days:02d}-{remaining_hours:02d}:{minutes:02d}"
    return(formatted_time)


def dorado_basecaller(toml_config):

    tool = "dorado_basecaller"
    output = toml_config["general"]["project_path"]
    email = toml_config["general"]["email"]
    genome = get_reference(toml_config["general"]["reference"])["fasta"]
    flowcells = toml_config["general"]["fc_dir_names"]
    
    # Iterate through each flowcell for basecalling and demultiplexing
    codes = []
    for flowcell in flowcells:
        reads = output + "/" + flowcell + "/reads/pod5"
        final = output + "/" + flowcell + "/alignments/"
        bam_dorado = final + flowcell + ".bam"

        # Simplify flowcell name with a shorter string
        code = flowcell.split('_')[-1]
        codes.append(code)

        cores = "8"
        memory = "64"

        # Get reads files size
        cmd = ["du", "-sh", "--apparent-size", "--block-size", "G", reads]
        result = subprocess.run(cmd, capture_output=True, text=True)

        size_str = result.stdout.split()[0].rstrip('G')

        # Scale required job time based on amount of data 
        hours = int(size_str) * 0.02
        formatted_time = format_time(hours)

        command = [TOOL_PATH + "main_pipelines/long-read/dorado-1.0.0-linux-x64/bin/dorado", 
                   "basecaller", "-v", "--device", "cuda:0", 
                   "--emit-moves", "--min-qscore", 
                   str(toml_config["dorado"]["min_q_score"]), 
                   "--reference", genome, 
                   "--sample-sheet", output + "/scripts/" + flowcell + ".csv", 
                   "--no-trim", "--kit-name", toml_config["general"]["kit"], 
                   "--mm2-opts", toml_config["dorado"]["mm2_opts"]
                   ]

        # Increases stringency if user defines it
        if toml_config["dorado"]["barcode_both_ends"] in ["true", "True", "yes", "Yes"]:
            command.extend(["--barcode-both-ends"])
        # Different model that includes base modification
        if "methylation" in toml_config["general"]["analysis"]:
            command.extend(["--modified-bases-models", TOOL_PATH + "main_pipelines/long-read/dorado_models/" + toml_config["dorado"]["modified_bases"], 
                            "--modified-bases-threshold", str(toml_config["dorado"]["modified_bases_threshold"])])
        # For transcriptomic data, when activated
        if "polya" in toml_config["general"]["analysis"]:
            command.extend(["--estimate-poly-a"])

        model = TOOL_PATH + "main_pipelines/long-read/dorado_models/" + toml_config['dorado']['model']
        command.extend([model, reads, ">", bam_dorado])

        command_str = " ".join(command)

        job = create_script(tool, cores, memory, formatted_time, output, email, command_str, flowcell)

        # Creates a variable job name for each flowcell (used for dependencies)
        var_name_bc = f"basecall_{code}"

        # Add slurm job to main.sh
        with open(output + "/scripts/main.sh", "a") as f:
            f.write("# Flowcell : " + flowcell + "\n")
            f.write(f"{var_name_bc}=$(sbatch --parsable {job})\n")

def dorado_demux(toml_config):

    tool = "dorado_demux"
    output = toml_config["general"]["project_path"]
    email = toml_config["general"]["email"]
    genome = get_reference(toml_config["general"]["reference"])["fasta"]
    flowcells = toml_config["general"]["fc_dir_names"]
    
    # Iterate through each flowcell for basecalling and demultiplexing
    codes = []
    for flowcell in flowcells:
        reads = output + "/" + flowcell + "/reads/pod5"
        final = output + "/" + flowcell + "/alignments/"
        bam_dorado = final + flowcell + ".bam"

        # Simplify flowcell name with a shorter string
        code = flowcell.split('_')[-1]
        codes.append(code)

        # Get reads files size
        cmd = ["du", "-sh", "--apparent-size", "--block-size", "G", reads]
        result = subprocess.run(cmd, capture_output=True, text=True)

        size_str = result.stdout.split()[0].rstrip('G')

        cores = "4"
        memory = "24"

        # Scale required job time based on amount of data
        hours = int(size_str) * 0.01
        formatted_time = format_time(hours)

        command = [TOOL_PATH + "main_pipelines/long-read/dorado-1.1.0-linux-x64/bin/dorado", 
                    "demux", "-vv", 
                    "--threads", cores, 
                    "--no-trim", 
                    "--output-dir", final, 
                    "--no-classify", bam_dorado, "\n\n"
                    ]
        command_str = " ".join(command)

        # Create slurm job
        job = create_script(tool, cores, memory, formatted_time, output, email, command_str, flowcell)

        # Different variable name for next set of dependencies
        var_name = f"demux_{code}"

        # Add slurm job to main.sh
        with open(output + "/scripts/main.sh", "a") as f:
            f.write(f"{var_name}=$(sbatch --parsable --dependency=afterok:${var_name_bc} {job})\n\n")

def samtools(toml_config):
    tool="samtools"
    cores="8"
    memory="32"
    time = "00-11:00"

    output = toml_config["general"]["project_path"]
    email = toml_config["general"]["email"]
    genome = get_reference(toml_config["general"]["reference"])["fasta"]
    flowcells = toml_config["general"]["fc_dir_names"]

    codes = []
    for flowcell in flowcells:
        code = flowcell.split('_')[-1]
        codes.append(code)
    dependencies = ":".join([f"$demux_{code}" for code in codes])
    command = ["python", "-u", 
                TOOL_PATH + "main_pipelines/long-read/LongReadSequencingONT/rename_bam.py", 
                toml_config["general"]["project_path"] + '/scripts/config_final.toml'
                ]
    command_str = " ".join(command)
    
    job = create_script(tool, cores, memory, time, output, email, command_str, "")

    # Add slurm job to main.sh
    with open(output + "/scripts/main.sh", "a") as f:
        f.write("# Rename, merge, sort and index bams")
        f.write(f"\nsamtools=$(sbatch --parsable --dependency=afterok:{dependencies} {job})\n")


def longreadsum(toml_config):
    tool="longreadsum"
    cores="4"
    memory="16"
    time = "00-23:00"

    output = toml_config["general"]["project_path"]
    email = toml_config["general"]["email"]
    genome = get_reference(toml_config["general"]["reference"])["fasta"]

    command_str = ""
    for name in toml_config["general"]["samples"]:
        command = ["apptainer", "run", 
                    TOOL_PATH + "main_pipelines/long-read/image_longreadsum.sif", 
                    "bam", 
                    "--threads", "8", 
                    "--ref", genome, 
                    "-Q", '"' + name + '_"', 
                    "-i", output + "/alignments/" + name + "_sorted.bam", 
                    "-o", output + "/qc"
                    ]

        if "methylation" in toml_config["general"]["analysis"]:
            command.extend(["--mod"])

        command.extend(["\n"])
        command.extend(["mv", output + "/qc/bam_summary.txt", output + "/qc/" + name + "_bam_summary.txt"])
        command.extend(["\n\n"])

        command_str += " ".join(command)

    job4 = create_script(tool, cores, memory, time, output, email, command_str, "")

    # Add slurm job to main.sh
    with open(output + "/scripts/main.sh", "a") as f:
        f.write("# QC\n")
        f.write(f"\nlongreadsum=$(sbatch --parsable --dependency=afterok:$samtools {job})\n")


def cleanup(toml_config):

    # Simple function to remove redundant files and cleanup structure 
    tool="cleanup"
    output = toml_config["general"]["project_path"]
    threads = "1"
    memory = "1"
    time = "00-00:30"
    email = toml_config["general"]["email"]
    flowcells = toml_config["general"]["fc_dir_names"]

    # Build cleanup commands
    commands = []

    # Move all logs to scripts/logs
    commands.append(f"mv {output}/*.log {output}/scripts/logs/")

    # Move all files in the main directory to scripts
    commands.append(f"mv -t {output}/scripts/ {output}/*.txt {output}/*.sh {output}/*.py")

    # Remove temp directories/files for each flowcell (MinKNOW-related OR redundant from merge)
    for flow in flowcells:
        commands.append(f"rm -rf {output}/{flow}/fastq_*")
        commands.append(f"rm -rf {output}/{flow}/bam_*")
        commands.append(f"rm -rf {output}/{flow}/alignments/*.bam")
        commands.append(f"rm -rf {output}/{flow}/alignments/*.bam.bai")
        commands.append(f"rm -f {output}/{flow}/main_reports/sequencing_summary*.txt")

    # Join all commands into a single string
    command_str = "\n".join(commands)

    # Create slurm job (Might be wasteful of server resources?)
    job = create_script(tool, threads, memory, time, output, email, command_str, "")

    done = []
    # For future - check if step is already done
    with open(output + "/scripts/steps_done.txt", "r") as f:
        for line in f:
            done.append(line.strip())

    # If QC is running, wait for it to finish
    if "qc" not in done:
        with open(output + "/scripts/main.sh", "a") as f:
            f.write("\n# Cleanup temporary files and logs\n")
            f.write(f"cleanup=$(sbatch --parsable --dependency=afterok:$longreadsum {job})\n")
    else:
        with open(output + "/scripts/main.sh", "a") as f:
            f.write("\n# Cleanup temporary files and logs\n")
            f.write(f"sbatch {job} \n")


# Launches main function
if __name__ == "__main__":
    main()
