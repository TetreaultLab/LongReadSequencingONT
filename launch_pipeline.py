import argparse
import toml
import time
import subprocess
import pandas as pd
import numpy as np
import os
from glob import glob
from datetime import date

# Global variables
TOOL_PATH = "/lustre09/project/6019267/shared/tools/"
DORADO = "main_pipelines/long-read/dorado-1.1.1-linux-x64/bin/dorado"


def main():
    parser = argparse.ArgumentParser(
        prog="PipelineLong",
        description="LongReadSequencing pipeline.",
    )

    parser.add_argument("config", type=str, help="Project config file, including path.")

    parser.add_argument(
        "--test",
        action="store_true",
        help="Enable testing: creates file without launching pipeline",
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
        if toml_config["general"]["kit"] != "SQK-LSK114":
            create_sample_sheet(toml_config)

        if not args.test:
            print(
                "\n\n\n!!! WARNING !!!\nIf you want to change the parameters: Press CTRL+C now!\nModify config_final.toml and launch_pipeline with that config file.\n\nOtherwise it will run with default parameters.\n\n"
            )
            time.sleep(10)

    else:
        toml_config = toml_config_initial

    # Get tool versions
    get_versions(output)

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
    if toml_config["general"]["kit"] == "SQK-LSK114":
        # Dorado Basecalling
        function_queue.append(dorado_samtools)

    else:
        # Dorado Basecalling
        function_queue.append(dorado_basecaller)

        # Dorado Basecalling
        function_queue.append(dorado_demux)

        # Renaming bams and running samtools merge, sort and index
        function_queue.append(samtools)

    # Call all other functions for downstream analysis
    # QC
    # function_queue.append(mosdepth)

    # EPI2ME
    function_queue.append(epi2me)

    # DNA specific
    if toml_config["general"]["seq_type"] == "WGS":
        # Repeat expansions
        if "repeats" in toml_config["general"]["analysis"]:
            function_queue.append(trgt)
            function_queue.append(strkit)

        # Methylation
        if "methylation" in toml_config["general"]["analysis"]:
            function_queue.append(ont_methyldmr_kit)

    # RNA specific
    if toml_config["general"]["seq_type"] == "RNA":
        # Splicing
        if "splicing" in toml_config["general"]["analysis"]:
            function_queue.append(flair)

        # Polyadenylation
        # if "polya" in toml_config["general"]["analysis"]:
        #     function_queue.append()

    # General variants analyses possible for DNA and RNA
    # SNPs
    if "SNP" in toml_config["general"]["analysis"]:
        function_queue.append(deepvariant)

    # # SVs
    if "SV" in toml_config["general"]["analysis"]:
        function_queue.append(cutesv)

    # # Annotation

    # # Phasing
    # if "phasing" in toml_config["general"]["analysis"]:
    #     function_queue.append()

    # Transfer back to projects
    function_queue.append(transfer)

    # Clean up
    function_queue.append(cleanup)

    # Create main.sh
    with open(output + "/scripts/main.sh", "w") as f:
        f.write("#!/bin/sh\n")

    # Calling main steps
    for func in function_queue:
        func(toml_config, done)

    # Check if in testing mode
    if args.test:
        print("\n\nTESTING MODE!\nThe pipeline will not be launched\n\n")
    else:
        # Call main.sh (Launch the pipeline)
        subprocess.run(["bash", output + "/scripts/main.sh"])


def create_config_final(filename):
    # Read initial config file
    with open(filename, "r") as f:
        toml_config = toml.load(f)

    # Set kit to variable
    kit = toml_config["general"]["kit"]

    # Check seq_type
    seq_type = toml_config["general"]["seq_type"]
    if (toml_config["general"]["seq_type"] == "Meth") or (
        "methylation" in toml_config["general"]["analysis"]
    ):
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
    toml_config["general"]["fc_dir_names"] = fc_dir_names

    # Making directory structure in scratch
    name = path.split("/")[-2].split("_", 1)[1]
    username = os.environ.get("USER")
    scratch = f"/lustre10/scratch/{username}/{name}/"

    directories = ["alignments", "results"] + fc_dir_names
    for d in directories:
        if not os.path.exists(scratch + "/" + d):
            os.makedirs(scratch + "/" + d)

    # Making directory structure in project
    directories = ["scripts", "scripts/logs", "alignments", "results", "qc"]
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
        reports = [
            "barcode_alignment",
            "final_summary",
            "pore_activity",
            "report",
            "sample_sheet",
            "sequencing_summary",
            "throughput",
            "output",
        ]
        for r in reports:
            f = glob(os.path.join(path + "/" + flow + "/", r + "*"))
            if f != []:  # if the file exist in the main directory
                for i in range(
                    0, len(f)
                ):  # if one or more file starts with the report name
                    f2 = f[i]
                    name = f2.split("/")[-1]
                    os.rename(
                        f2, os.path.join(path + "/" + flow + "/main_reports", name)
                    )  # move to main_reports

        # Move pod5 to reads directory
        if not os.path.exists(path + "/" + flow + "/reads/pod5"):
            os.rename(
                os.path.join(path, flow, "pod5"),
                os.path.join(path, flow, "reads", "pod5"),
            )
        if os.path.exists(path + "/" + flow + "/pod5_skip"):
            for file in os.listdir(os.path.join(path, flow, "pod5_skip")):
                os.rename(
                    os.path.join(path, flow, "pod5_skip", file),
                    os.path.join(path, flow, "reads", "pod5", file),
                )

    # Add Dorado options
    ## general options
    toml_config["dorado"] = {}
    toml_config["dorado"]["min_q_score"] = 9
    toml_config["dorado"]["barcode_both_ends"] = "False"
    toml_config["dorado"]["model"] = "dna_r10.4.1_e8.2_400bps_sup@v5.2.0"

    ## parameters depending on sequencing type
    if seq_type == "WGS":
        toml_config["dorado"]["mm2_opts"] = '"-x lr:hq"'
        if kit not in ["SQK-RBK114-24", "SQK-NBD114-24", "SQK-LSK114"]:
            raise Exception(
                "Error: Wrong Kit for WGS. Options are SQK-RBK114-24, SQK-NBD114-24, SQK-LSK114"
            )

    if seq_type == "RNA":
        toml_config["dorado"]["mm2_opts"] = '"-x splice:hq"'
        if kit not in ["SQK-PCB114-24"]:
            raise Exception(
                "Error: Wrong Kit for Whole Transcriptome. Options are SQK-PCB114-24"
            )

    if methylation_status:
        toml_config["dorado"]["modified_bases"] = (
            "dna_r10.4.1_e8.2_400bps_sup@v5.2.0_5mC_5hmC@v2"
        )
        toml_config["dorado"]["modified_bases_threshold"] = 0.0

    if seq_type == "Targeted":
        toml_config["dorado"]["mm2_opts"] = '"-x splice --junc-bed anno.bed12"'
        if kit not in ["SQK-NBD114-24"]:
            raise Exception(
                "Error: Wrong Kit for Targeted Sequencing. Options are SQK-NBD114-24"
            )

    # Add mosdepth options
    ## general options
    toml_config["mosdepth"] = {}
    toml_config["mosdepth"]["bins"] = 25000  # Should re-use with wf-human-variation
    toml_config["mosdepth"]["thresholds"] = "1,5,10,20,30"

    # Add TRGT options
    toml_config["trgt"] = {}
    toml_config["trgt"]["gene_interest"] = ""
    toml_config["trgt"]["motif"] = ""

    # Add next tool options

    # Save new config file
    with open(
        toml_config["general"]["project_path"] + "/scripts/config_final.toml", "w"
    ) as f:
        toml.dump(toml_config, f)

    return toml_config


def create_sample_sheet(toml_config):
    path = toml_config["general"]["project_path"].rstrip("/")
    parts = path.split(os.sep)

    if len(parts) < 2:
        raise ValueError(
            f"Project path '{path}' is too short to parse project/date info"
        )

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
            d = {
                "flow_cell_id": flow_cell_id,
                "experiment_id": project_name,
                "kit": kit,
                "alias": samples,
                "type": conditions,
                "barcode": barcode,
            }
        else:
            d = {
                "flow_cell_id": flow_cell_id,
                "experiment_id": project_name,
                "kit": kit,
                "alias": samples,
                "type": conditions,
                "barcode": range(barcode, barcode + len(samples)),
            }

        # Create the sample sheet and save to .csv
        df = pd.DataFrame(data=d)
        df["barcode"] = "barcode" + df["barcode"].astype(int).astype(str).str.zfill(2)

        df.to_csv(path + "/scripts/" + flowcell + ".csv", sep=",", index=False)


def get_versions(output):
    today = date.today()
    # subprocess.run(
    #     ["module", "load", "samtools", "apptainer/1.4.5", "nextflow/25.04.6"]
    # )

    cutesv = subprocess.run(
        [
            "apptainer",
            "run",
            TOOL_PATH + "variants/SVs/cutesv/cutesv.sif",
            "cuteSV",
            "--version",
        ],
        capture_output=True,
        text=True,
    )

    deepvariant = subprocess.run(
        [
            "apptainer",
            "run",
            TOOL_PATH + "variants/SNPs/deepvariant/deepvariant.sif",
            "/opt/deepvariant/bin/run_deepvariant",
            "--version",
        ],
        capture_output=True,
        text=True,
    )

    dorado = subprocess.run(
        [TOOL_PATH + DORADO, "--version"], capture_output=True, text=True
    )

    epi2me = subprocess.run(
        [
            "nextflow",
            "run",
            "epi2me-labs/wf-human-variation",
            "--version",
        ],
        capture_output=True,
        text=True,
    )

    flair = subprocess.run(
        [
            "apptainer",
            "run",
            TOOL_PATH + "splicing/flair/flair.sif",
            "flair",
            "--version",
        ],
        capture_output=True,
        text=True,
    )

    intron_prospector = subprocess.run(
        [
            "apptainer",
            "run",
            TOOL_PATH + "splicing/flair/flair.sif",
            TOOL_PATH + "splicing/flair/intronProspector/bin/intronProspector",
            "--version",
        ],
        capture_output=True,
        text=True,
    )

    ont_methyldmr_kit = "ont-methylDMR-kit version 3.2"

    samtools = subprocess.run(["samtools", "--version"], capture_output=True, text=True)

    strkit = subprocess.run(
        [
            "apptainer",
            "run",
            TOOL_PATH + "repeat_expansions/STRkit/strkit.sif",
            "strkit",
            "--version",
        ],
        capture_output=True,
        text=True,
    )

    trgt = subprocess.run(
        [
            "apptainer",
            "run",
            TOOL_PATH + "repeat_expansions/TRGT/trgt.sif",
            "trgt",
            "--version",
        ],
        capture_output=True,
        text=True,
    )

    with open(f"{output}/scripts/logs/versions_{today}.txt", "w") as f:
        f.write(f"LRS pipeline tool versions {today}")
        f.write(cutesv)
        f.write(deepvariant)
        f.write(dorado)
        f.write(epi2me)
        f.write(flair)
        f.write(intron_prospector)
        f.write(ont_methyldmr_kit)
        f.write(samtools)
        f.write(strkit)
        f.write(trgt)


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


def create_script(tool, cores, memory, time, output, email, command, flowcell):
    steps_done = output + "/scripts/steps_done.txt"
    name = output.rstrip("/").split("/")[-2].split("_", 1)[1]

    # Enables creating a script per flowcell, or a single script if "" is added as the argument
    if flowcell != "":
        code = flowcell.split("_")[-1]
        job = output + "/scripts/" + tool + "_" + flowcell + ".slurm"

        # Uses a slurm template for each job script
        with open(
            TOOL_PATH
            + "main_pipelines/long-read/LongReadSequencingONT/template_sbatch.txt",
            "r",
        ) as f:
            slurm = f.read()
            # Because dorado utilizes a GPU, it needs different options
            if tool == "dorado_basecaller":
                slurm_filled = slurm.format(
                    cores,
                    "#SBATCH --gpus=h100:1",
                    memory,
                    time,
                    tool,
                    flowcell,
                    "log",
                    "log",
                    "def",
                    email,
                    name,
                )

            elif tool == "dorado_demux":
                slurm_filled = slurm.format(
                    cores,
                    "",
                    memory,
                    time,
                    tool,
                    flowcell,
                    "out",
                    "log",
                    "rrg",
                    email,
                    name,
                )

            # Most tools will be CPU-dependent
            else:
                slurm_filled = slurm.format(
                    cores,
                    "",
                    memory,
                    time,
                    tool,
                    flowcell,
                    "log",
                    "log",
                    "rrg",
                    email,
                    name,
                )
                slurm_filled += "module load StdEnv/2023 apptainer samtools\n"
                slurm_filled += (
                    "source "
                    + TOOL_PATH
                    + "main_pipelines/long-read/launch_pipeline_env/bin/activate"
                )

            slurm_filled += "\n#\n### Calling " + tool + " - " + flowcell + "\n#\n"
            slurm_filled += command
            slurm_filled += "\n"

            # Keep track of completed steps
            slurm_filled += (
                f'if [ $? -eq 0 ]; then echo "{tool}_{code}" >> "{steps_done}"; fi\n\n'
            )

    # This is for tools running on all the data at once
    else:
        if tool == "samtools":
            job = output + "/scripts/" + tool + ".slurm"

            with open(
                TOOL_PATH
                + "main_pipelines/long-read/LongReadSequencingONT/template_sbatch.txt",
                "r",
            ) as f:
                slurm = f.read()
                # For all other tools
                slurm_filled = slurm.format(
                    cores,
                    "#SBATCH --ntasks=5",
                    memory,
                    time,
                    tool,
                    "run",
                    "log",
                    "log",
                    "rrg",
                    email,
                    name,
                )

                # Add enviroment loading commands
                slurm_filled += "module load StdEnv/2023 apptainer samtools\n"
                slurm_filled += (
                    "source "
                    + TOOL_PATH
                    + "main_pipelines/long-read/launch_pipeline_env/bin/activate\n"
                )

                slurm_filled += "\n#\n### Calling " + tool + "\n#\n"
                slurm_filled += command
                slurm_filled += "\n"

                # Keep track of completed steps
                slurm_filled += (
                    f'if [ $? -eq 0 ]; then echo "{tool}" >> "{steps_done}"; fi\n\n'
                )
        else:
            job = output + "/scripts/" + tool + ".slurm"

            with open(
                TOOL_PATH
                + "main_pipelines/long-read/LongReadSequencingONT/template_sbatch.txt",
                "r",
            ) as f:
                slurm = f.read()
                # For all other tools
                slurm_filled = slurm.format(
                    cores,
                    "",
                    memory,
                    time,
                    tool,
                    "run",
                    "log",
                    "log",
                    "rrg",
                    email,
                    name,
                )

                # Add enviroment loading commands
                slurm_filled += "module load StdEnv/2023 apptainer samtools\n"
                slurm_filled += (
                    "source "
                    + TOOL_PATH
                    + "main_pipelines/long-read/launch_pipeline_env/bin/activate\n"
                )

                slurm_filled += "\n#\n### Calling " + tool + "\n#\n"
                slurm_filled += command
                slurm_filled += "\n"

                # Keep track of completed steps
                slurm_filled += (
                    f'if [ $? -eq 0 ]; then echo "{tool}" >> "{steps_done}"; fi\n\n'
                )

    with open(job, "w") as o:
        o.write(slurm_filled)

        return job


def format_time(hours):
    # Optimization of time allocation (buffers within AllianceCan priority)
    if hours < 2.5:  # Rule 1: For hours < 2.5
        if hours * 2 > 3:  # Only weakness I would see is 2.0-2.5h (may timeout)
            hours = 3.0
        else:
            hours = np.ceil(hours * 2)
    elif 2.5 <= hours < 6:  # Rule 2: For 2.5 ≤ hours < 10
        hours = np.ceil(hours + 1)  # +1 hour buffer
    elif 6 <= hours < 10:
        hours = np.ceil(hours + 2)  # +2 hours buffer
    elif 10 <= hours < 15:  # Rule 3: For 10 ≤ hours < 20
        hours = np.ceil(hours + 3)  # +3 hours buffer
    elif 15 <= hours < 20:
        hours = np.ceil(hours + 4)  # +4 hours buffer
    elif 20 <= hours < 66:  # Rule 4: For 20 ≤ hours < 66
        hours = np.ceil(hours + 6)  # +6 hours buffer
    elif 66 <= hours < 158:  # Rule 5: For 66 ≤ hours < 158
        hours = np.ceil(hours + 10)  # +10 hours buffer
    elif hours >= 158:  # Rule 6: For 158 ≤ hours
        hours = 168  # Max SLURM time: 7 days
        print(
            "!!! Warning: The amount of data suggests the job may time out after 7 days."
        )

    # Formatting variables
    days = int(hours // 24)
    remaining_hours = int(hours % 24)
    minutes = int((hours % 1) * 60)

    # Format as DD-HH:MM
    formatted_time = f"{days:02d}-{remaining_hours:02d}:{minutes:02d}"
    return formatted_time


def dorado_basecaller(toml_config, done):
    tool = "dorado_basecaller"
    output = toml_config["general"]["project_path"]
    email = toml_config["general"]["email"]
    genome = get_reference(toml_config["general"]["reference"])["fasta"]
    flowcells = toml_config["general"]["fc_dir_names"]
    name = output.rstrip("/").split("/")[-2].split("_", 1)[1]
    username = os.environ.get("USER")

    # Iterate through each flowcell for basecalling
    for flowcell in flowcells:
        reads = output + "/" + flowcell + "/reads/pod5"
        final = f"/lustre10/scratch/{username}/{name}/{flowcell}/"
        bam_dorado = final + flowcell + ".bam"

        cores = "8"
        memory = "64"

        # Get reads files size
        cmd = ["du", "-sh", "--apparent-size", "--block-size", "G", reads]
        result = subprocess.run(cmd, capture_output=True, text=True)

        size_str = result.stdout.split()[0].rstrip("G")

        code = flowcell.split("_")[-1]

        # Scale required job time based on amount of data
        hours = int(size_str) * 0.04
        formatted_time = format_time(hours)

        command = [
            TOOL_PATH + DORADO,
            "basecaller",
            "-v",
            "--device",
            "cuda:0",
            "--min-qscore",
            str(toml_config["dorado"]["min_q_score"]),
            "--reference",
            genome,
            "--sample-sheet",
            output + "/scripts/" + flowcell + ".csv",
            "--no-trim",
            "--kit-name",
            toml_config["general"]["kit"],
            "--mm2-opts",
            toml_config["dorado"]["mm2_opts"],
        ]

        # Increases stringency if user defines it
        if toml_config["dorado"]["barcode_both_ends"] in ["true", "True", "yes", "Yes"]:
            command.extend(["--barcode-both-ends"])
        # Different model that includes base modification
        if "methylation" in toml_config["general"]["analysis"]:
            command.extend(
                [
                    "--modified-bases-models",
                    TOOL_PATH
                    + "main_pipelines/long-read/dorado_models/"
                    + toml_config["dorado"]["modified_bases"],
                    "--modified-bases-threshold",
                    str(toml_config["dorado"]["modified_bases_threshold"]),
                ]
            )
        # For transcriptomic data, when activated
        if "polya" in toml_config["general"]["analysis"]:
            command.extend(["--estimate-poly-a"])

        model = (
            TOOL_PATH
            + "main_pipelines/long-read/dorado_models/"
            + toml_config["dorado"]["model"]
        )
        command.extend([model, reads, ">", bam_dorado])

        command_str = " ".join(command)

        job = create_script(
            tool, cores, memory, formatted_time, output, email, command_str, flowcell
        )

        # Creates a variable job name for each flowcell (used for dependencies)
        var_name_bc = f"dorado_basecaller_{code}"

        # Add slurm job to main.sh
        if var_name_bc not in done:
            print("To-Do: " + var_name_bc)
            with open(output + "/scripts/main.sh", "a") as f:
                f.write(f"\n# Dorado Basecall for flowcell : {flowcell}")
                f.write(f"\n{var_name_bc}=$(sbatch --parsable {job})\n")
        else:
            print("Done: " + var_name_bc)


def dorado_demux(toml_config, done):
    tool = "dorado_demux"
    output = toml_config["general"]["project_path"]
    email = toml_config["general"]["email"]
    flowcells = toml_config["general"]["fc_dir_names"]
    name = output.rstrip("/").split("/")[-2].split("_", 1)[1]
    username = os.environ.get("USER")

    cores = "2"
    memory = "8"

    # Iterate through each flowcell for demultiplexing
    for flowcell in flowcells:
        reads = output + "/" + flowcell + "/reads/pod5"
        final = f"/lustre10/scratch/{username}/{name}/{flowcell}/"
        bam_dorado = final + flowcell + ".bam"

        # Get reads files size
        cmd = ["du", "-sh", "--apparent-size", "--block-size", "G", reads]
        result = subprocess.run(cmd, capture_output=True, text=True)

        size_str = result.stdout.split()[0].rstrip("G")

        # Scale required job time based on amount of data
        hours = int(size_str) * 0.02
        formatted_time = format_time(hours)

        code = flowcell.split("_")[-1]

        command = [
            TOOL_PATH + DORADO,
            "demux",
            "-vv",
            "--threads",
            cores,
            "--no-trim",
            "--output-dir",
            final,
            "--no-classify",
            bam_dorado,
            "\n\n",
        ]
        command_str = " ".join(command)

        # Create slurm job
        job = create_script(
            tool, cores, memory, formatted_time, output, email, command_str, flowcell
        )

        # Different variable name for next set of dependencies
        var_name = f"dorado_demux_{code}"
        var_name_bc = f"dorado_basecaller_{code}"

        # Add slurm job to main.sh
        if var_name not in done:  # demux not done
            print("To-Do: " + var_name)
            if var_name_bc not in done:  # basecall and demux not done
                with open(output + "/scripts/main.sh", "a") as f:
                    f.write(f"\n# Dorado Demux for flowcell : {flowcell}")
                    f.write(
                        f"\n{var_name}=$(sbatch --parsable --dependency=afterok:${var_name_bc} {job})\n"
                    )
            else:  # demux not done but basecall done
                with open(output + "/scripts/main.sh", "a") as f:
                    f.write(f"\n# Dorado Demux for flowcell : {flowcell}")
                    f.write(f"\n{var_name}=$(sbatch --parsable {job})\n")
        else:  # demux done
            if (
                var_name_bc not in done
            ):  # demux done but not basecall, redo demux after basecall
                print("To-Do: " + var_name)
                with open(output + "/scripts/main.sh", "a") as f:
                    f.write(f"\n# Dorado Demux for flowcell : {flowcell}")
                    f.write(
                        f"\n{var_name}=$(sbatch --parsable --dependency=afterok:${var_name_bc} {job})\n"
                    )
            else:  # demux done and basecall done, sucess
                print("Done: " + var_name)


def dorado_samtools(toml_config, done):
    tool = "dorado_basecaller"

    output = toml_config["general"]["project_path"]
    email = toml_config["general"]["email"]
    genome = get_reference(toml_config["general"]["reference"])["fasta"]
    flowcell = toml_config["general"]["fc_dir_names"][0]
    code = flowcell.split("_")[-1]
    sample = toml_config["general"]["samples"][0]
    name = output.rstrip("/").split("/")[-2].split("_", 1)[1]
    username = os.environ.get("USER")
    tmpdir = os.environ.get("SLURM_TMPDIR")
    sorted_bam = f"$SLURM_TMPDIR/{sample}_sorted.bam"

    reads = output + "/" + flowcell + "/reads/pod5"
    final = f"/lustre10/scratch/{username}/{name}/alignments/"
    bam_dorado = f"{final}{sample}.bam"

    cores = "8"
    memory = "64"

    # Get reads files size
    cmd = ["du", "-sh", "--apparent-size", "--block-size", "G", reads]
    result = subprocess.run(cmd, capture_output=True, text=True)

    size_str = result.stdout.split()[0].rstrip("G")

    # Scale required job time based on amount of data
    hours = int(size_str) * 0.04
    formatted_time = format_time(hours)

    command = [
        TOOL_PATH + DORADO,
        "basecaller",
        "-v",
        "--device",
        "cuda:0",
        "--min-qscore",
        str(toml_config["dorado"]["min_q_score"]),
        "--reference",
        genome,
        "--no-trim",
        "--mm2-opts",
        toml_config["dorado"]["mm2_opts"],
    ]

    # Increases stringency if user defines it
    if toml_config["dorado"]["barcode_both_ends"] in ["true", "True", "yes", "Yes"]:
        command.extend(["--barcode-both-ends"])
    # Different model that includes base modification
    if "methylation" in toml_config["general"]["analysis"]:
        command.extend(
            [
                "--modified-bases-models",
                TOOL_PATH
                + "main_pipelines/long-read/dorado_models/"
                + toml_config["dorado"]["modified_bases"],
                "--modified-bases-threshold",
                str(toml_config["dorado"]["modified_bases_threshold"]),
            ]
        )
    # For transcriptomic data, when activated
    if "polya" in toml_config["general"]["analysis"]:
        command.extend(["--estimate-poly-a"])

    model = (
        TOOL_PATH
        + "main_pipelines/long-read/dorado_models/"
        + toml_config["dorado"]["model"]
    )
    command.extend([model, reads, ">", bam_dorado])

    # Add samtools sort and index to command
    command.extend(
        [
            "\n\nml samtools",
            "\n\n#Sort\n",
            "samtools",
            "sort",
            "--threads",
            "7",
            "-m",
            "4G",
            "-o",
            sorted_bam,
            bam_dorado,
        ]
    )

    # transfer
    command.extend(["\n\n#Transfer\n", "cp", sorted_bam, f"{final}{sample}_sorted.bam"])

    # index
    command.extend(
        [
            "\n\n#Index\n",
            "samtools",
            "index",
            "--threads",
            "7",
            f"{final}{sample}_sorted.bam",
            "\n",
        ]
    )

    command_str = " ".join(command)

    job = create_script(
        tool, cores, memory, formatted_time, output, email, command_str, flowcell
    )

    # Add slurm job to main.sh
    if f"{tool}_{code}" not in done:
        print("To-Do: dorado_samtools")
        with open(output + "/scripts/main.sh", "a") as f:
            f.write(f"\n# Dorado Basecall and Samtools for sample : {sample}")
            f.write(f"\nsamtools_{sample}=$(sbatch --parsable {job})\n")
    else:
        print("Done: dorado_samtools")


def samtools(toml_config, done):
    tool = "samtools"

    output = toml_config["general"]["project_path"]
    email = toml_config["general"]["email"]
    flowcells = toml_config["general"]["fc_dir_names"]
    samples = toml_config["general"]["samples"]
    config = toml_config["general"]["project_path"] + "/scripts/config_final.toml"

    codes = []
    for flowcell in flowcells:
        code = flowcell.split("_")[-1]
        codes.append(code)

    all_fc = [f"dorado_demux_{code}" for code in codes]
    done_fc = [x for x in done if x.startswith("dorado_demux")]
    to_dos = [x for x in all_fc if x not in done_fc]

    for sample in samples:
        job = output + "/scripts/" + tool + "_" + sample + ".slurm"
        with open(
            TOOL_PATH
            + "main_pipelines/long-read/LongReadSequencingONT/template_samtools.txt",
            "r",
        ) as f:
            slurm = f.read()
            slurm_filled = slurm.format(sample, email, output, config, sample)

            with open(job, "w") as o:
                o.write(slurm_filled)

        samtools_name = f"samtools_{sample}"

        if len(to_dos) > 0:
            print("To-Do: " + samtools_name)
            dependencies = ":".join([f"${code}" for code in to_dos])
            with open(output + "/scripts/main.sh", "a") as f:
                f.write("\n# Rename, merge, sort and index bams")
                f.write(
                    f"\n{samtools_name}=$(sbatch --parsable --dependency=afterok:{dependencies} {job})\n"
                )

            # remove samtools from done so all subsequent jobs run after samtools
            if samtools_name in done:
                done.remove(samtools_name)
        else:
            # If all dorado_demux are done but samtools has not run yet
            if samtools_name not in done:
                print("To-Do: " + samtools_name)
                with open(output + "/scripts/main.sh", "a") as f:
                    f.write("\n# Rename, merge, sort and index bams for {sample}")
                    f.write(f"\n{samtools_name}=$(sbatch --parsable {job})\n")
            else:
                # All dorado_demux are done and samtools is done
                print("Done: " + samtools_name)


def longReadSum(toml_config, done):
    tool = "longReadSum"
    cores = "4"
    memory = "16"
    time = "00-23:00"

    output = toml_config["general"]["project_path"]
    email = toml_config["general"]["email"]
    genome = get_reference(toml_config["general"]["reference"])["fasta"]

    command_str = ""
    for name in toml_config["general"]["samples"]:
        command = [
            "apptainer",
            "run",
            TOOL_PATH + "main_pipelines/long-read/image_longreadsum.sif",
            "bam",
            "--threads",
            cores,
            "--ref",
            genome,
            "-Q",
            '"' + name + '_"',
            "-i",
            output + "/alignments/" + name + "_sorted.bam",
            "-o",
            output + "/qc",
        ]

        if "methylation" in toml_config["general"]["analysis"]:
            command.extend(["--mod"])

        command.extend(["\n"])
        command.extend(
            [
                "mv",
                output + "/qc/bam_summary.txt",
                output + "/qc/" + name + "_bam_summary.txt",
            ]
        )
        command.extend(["\n\n"])

        command_str += " ".join(command)

    job = create_script(tool, cores, memory, time, output, email, command_str, "")

    # Add slurm job to main.sh
    if "samtools" not in done:
        print("To-Do: " + tool)
        with open(output + "/scripts/main.sh", "a") as f:
            f.write("\n# QC")
            f.write(
                f"\nlongreadsum=$(sbatch --parsable --dependency=afterok:$samtools {job})\n"
            )
    else:
        if tool not in done:
            print("To-Do: " + tool)
            with open(output + "/scripts/main.sh", "a") as f:
                f.write("\n# QC")
                f.write(f"\nlongreadsum=$(sbatch --parsable {job})\n")
        else:
            print("Done: " + tool)


def mosdepth(toml_config, done):
    # As a module until nextflow is usable
    tool = "mosdepth"
    output = toml_config["general"]["project_path"]
    flowcells = toml_config["general"]["fc_dir_names"]
    threads = "4"
    memory = "8"
    email = toml_config["general"]["email"]

    dirs = [f"{output}/{fc}/reads/pod5" for fc in flowcells]
    cmd = ["du", "-sh", "--apparent-size", "--block-size", "G", "--total"] + dirs
    result = subprocess.run(cmd, capture_output=True, text=True)

    for line in result.stdout.splitlines():
        if line.endswith("total"):
            size = line.split()[0]

    size_str = size.rstrip("G")
    hours = int(size_str) * 0.003
    formatted_time = format_time(hours)

    command_str = ""
    for name in toml_config["general"]["samples"]:
        input_file = output + "/alignments/" + name + "_sorted.bam"
        # Main mosdepth function
        command = [
            "apptainer",
            "run",
            TOOL_PATH + "others/mosdepth/mosdepth.sif",
            "mosdepth",
            "-t",
            threads,
            "-n",
            "-x",
            "-b",
            str(toml_config["mosdepth"]["bins"]),
            "-T",
            str(toml_config["mosdepth"]["thresholds"]),
            output + "/qc/" + name,
            input_file,
        ]
        command_str += " ".join(command) + "\n"
        # Added visualization function
        command2 = [
            "python",
            "-u",
            TOOL_PATH + "others/mosdepth/SummarizedMosdepth.py",
            "-p",
            output + "/qc/" + name,
            "--bins",
            str(toml_config["mosdepth"]["bins"]),
            "--thresholds",
            str(toml_config["mosdepth"]["thresholds"]),
        ]
        command_str += " ".join(command2) + "\n"

    # Make a report with all the plots generated (Only once per run)
    command3 = [
        "python",
        "-u",
        TOOL_PATH + "others/mosdepth/mosdepth_report.py",
        "-i",
        output + "/qc",
    ]
    command_str += " ".join(command3) + "\n"

    job = create_script(
        tool, threads, memory, formatted_time, output, email, command_str, ""
    )

    # Add slurm job to main.sh
    if "samtools" not in done:
        print("To-Do: " + tool)
        with open(output + "/scripts/main.sh", "a") as f:
            f.write("\n# Mosdepth")
            f.write(
                f"\nmosdepth=$(sbatch --parsable --dependency=afterok:$samtools {job})\n"
            )

    else:
        if tool not in done:
            print("To-Do: " + tool)
            with open(output + "/scripts/main.sh", "a") as f:
                f.write("\n# Mosdepth")
                f.write(f"\nmosdepth=$(sbatch --parsable {job})\n")
        else:
            print("Done: " + tool)


def toulligqc(toml_config, done):
    tool = "toulligqc"
    output = toml_config["general"]["project_path"]
    flowcells = toml_config["general"]["fc_dir_names"]
    threads = "4"
    memory = "8"
    email = toml_config["general"]["email"]

    dirs = [f"{output}/{fc}/reads/pod5" for fc in flowcells]
    cmd = ["du", "-sh", "--apparent-size", "--block-size", "G", "--total"] + dirs
    result = subprocess.run(cmd, capture_output=True, text=True)

    for line in result.stdout.splitlines():
        if line.endswith("total"):
            size = line.split()[0]

    size_str = size.rstrip("G")
    hours = int(size_str) * 0.003
    formatted_time = format_time(hours)

    subprocess.run(["rm", "-r", output + "/qc/"])

    seq_summary = []
    for fc in flowcells:
        pattern = f"{output}/{fc}/main_reports/sequencing_summary_*.txt"
        for f in glob(pattern):
            seq_summary.append(f"-a {f}")

    samplesheets = []
    for fc in flowcells:
        pattern = f"{output}/{fc}/main_reports/sample_sheet_*.csv"
        for f in glob(pattern):
            samplesheets.append(f"-s {f}")

    bams = []
    for name in toml_config["general"]["samples"]:
        input_file = output + "/alignments/" + name + "_sorted.bam"
        bams.append(f"-u {input_file}")

    command = (
        [
            "apptainer",
            "run",
            TOOL_PATH + "main_pipelines/long-read/toulligqc.sif",
            "toulligqc",
            "--output-directory",
            output + "/qc/",
            "--report-name",
            "toulligqc",
            "--thread",
            threads,
            "--barcoding",
        ]
        + seq_summary
        # + samplesheets
        + bams
    )

    command_str = " ".join(command)
    print(command_str)

    job = create_script(
        tool, threads, memory, formatted_time, output, email, command_str, ""
    )

    if "samtools" not in done:
        print("To-Do: " + tool)
        with open(output + "/scripts/main.sh", "a") as f:
            f.write("\n# ToulligQC")
            f.write(
                f"\ntoulligqc=$(sbatch --parsable --dependency=afterok:$samtools {job})\n"
            )

    else:
        if tool not in done:
            print("To-Do: " + tool)
            with open(output + "/scripts/main.sh", "a") as f:
                f.write("\n# ToulligQC")
                f.write(f"\ntoulligqc=$(sbatch --parsable {job})\n")
        else:
            print("Done: " + tool)


def epi2me(toml_config, done):
    tool = "epi2me"
    cores = "24"
    memory = "96"
    time = "00-01:00"

    output = toml_config["general"]["project_path"]
    email = toml_config["general"]["email"]
    model = toml_config["dorado"]["model"]

    genome = get_reference(toml_config["general"]["reference"])["fasta"]
    analysis = toml_config["general"]["analysis"]
    name = output.rstrip("/").split("/")[-2].split("_", 1)[1]

    todo = ""
    arg_map = {
        "methylation": "--mod",
        "SNP": "--snp",
        "SV": "--sv",
        "CNV": "--cnv",
        "repeats": "--str",
        "phasing": "--phased",
    }

    for item in analysis:
        if item in arg_map:
            todo += " " + arg_map[item]

    for sample in toml_config["general"]["samples"]:
        job = output + "/scripts/" + tool + "_" + sample + ".slurm"
        with open(
            TOOL_PATH
            + "main_pipelines/long-read/LongReadSequencingONT/template_epi2me.txt",
            "r",
        ) as f:
            slurm = f.read()
            slurm_filled = slurm.format(
                cores,
                memory,
                time,
                tool,
                sample,
                email,
                name,
                output,
                todo,
                model,
                genome,
            )

            with open(job, "w") as o:
                o.write(slurm_filled)

        samtools_name = f"samtools_{sample}"
        epi_name = f"epi2me_{sample}"

        if samtools_name not in done:
            print("To-Do: " + epi_name)
            with open(output + "/scripts/main.sh", "a") as f:
                f.write(f"\n# Epi2me workflow human variation for {sample}")
                f.write(
                    f"\n{epi_name}=$(sbatch --parsable --dependency=afterok:${samtools_name} {job})\n"
                )
        else:
            if epi_name not in done:
                print("To-Do: " + epi_name)
                with open(output + "/scripts/main.sh", "a") as f:
                    f.write(f"\n# Epi2me workflow human variation for {sample}")
                    f.write(f"\n{epi_name}=$(sbatch --parsable {job})\n")
            else:
                print("Done: " + epi_name)


def trgt(toml_config, done):
    tool = "trgt"

    output = toml_config["general"]["project_path"]
    email = toml_config["general"]["email"]
    genome = get_reference(toml_config["general"]["reference"])["fasta"]
    name = output.rstrip("/").split("/")[-2].split("_", 1)[1]
    samples = toml_config["general"]["samples"]
    str_samples = " ".join(samples)

    job = output + "/scripts/" + tool + ".slurm"
    with open(
        TOOL_PATH + "main_pipelines/long-read/LongReadSequencingONT/template_trgt.txt",
        "r",
    ) as f:
        slurm = f.read()
        slurm_filled = slurm.format(email, output, genome, str_samples, name)

        with open(job, "w") as o:
            o.write(slurm_filled)

    all_fc = [f"samtools_{s}" for s in samples]
    done_fc = [x for x in done if x.startswith("samtools")]
    to_dos = [x for x in all_fc if x not in done_fc]

    dependencies = ":".join([f"${code}" for code in to_dos])

    if len(to_dos) > 0:
        print("To-Do: " + tool)
        with open(output + "/scripts/main.sh", "a") as f:
            f.write("\n# TRGT")
            f.write(
                f"\ntrgt=$(sbatch --parsable --dependency=afterok:{dependencies} {job})\n"
            )
    else:
        if tool not in done:
            print("To-Do: " + tool)
            with open(output + "/scripts/main.sh", "a") as f:
                f.write("\n# TRGT")
                f.write(f"\ntrgt=$(sbatch --parsable {job})\n")
        else:
            print("Done: " + tool)


def strkit(toml_config, done):
    tool = "strkit"

    output = toml_config["general"]["project_path"]
    email = toml_config["general"]["email"]
    genome = get_reference(toml_config["general"]["reference"])["fasta"]
    name = output.rstrip("/").split("/")[-2].split("_", 1)[1]
    samples = toml_config["general"]["samples"]
    str_samples = " ".join(samples)

    job = output + "/scripts/" + tool + ".slurm"
    with open(
        TOOL_PATH
        + "main_pipelines/long-read/LongReadSequencingONT/template_strkit.txt",
        "r",
    ) as f:
        slurm = f.read()
        slurm_filled = slurm.format(email, output, genome, str_samples, name)

        with open(job, "w") as o:
            o.write(slurm_filled)

    all_fc = [f"samtools_{s}" for s in samples]
    done_fc = [x for x in done if x.startswith("samtools")]
    to_dos = [x for x in all_fc if x not in done_fc]

    dependencies = ":".join([f"${code}" for code in to_dos])

    if len(to_dos) > 0:
        print("To-Do: " + tool)
        with open(output + "/scripts/main.sh", "a") as f:
            f.write("\n# STRkit")
            f.write(
                f"\nstrkit=$(sbatch --parsable --dependency=afterok:{dependencies} {job})\n"
            )
    else:
        if tool not in done:
            print("To-Do: " + tool)
            with open(output + "/scripts/main.sh", "a") as f:
                f.write("\n# STRkit")
                f.write(f"\nstrkit=$(sbatch --parsable {job})\n")
        else:
            print("Done: " + tool)


def ont_methyldmr_kit(toml_config, done):
    tool = "ont_methyldmr_kit"

    output = toml_config["general"]["project_path"]
    email = toml_config["general"]["email"]
    name = output.rstrip("/").split("/")[-2].split("_", 1)[1]
    username = os.environ.get("USER")
    res = f"/lustre10/scratch/{username}/{name}/results"
    samples = toml_config["general"]["samples"]
    str_samples = " ".join(samples)

    job = output + "/scripts/" + tool + ".slurm"
    with open(
        TOOL_PATH
        + "main_pipelines/long-read/LongReadSequencingONT/template_ont-methyldmr-kit.txt",
        "r",
    ) as f:
        slurm = f.read()
        slurm_filled = slurm.format(email, output, str_samples, res)

        with open(job, "w") as o:
            o.write(slurm_filled)

    all_fc = [f"epi2me_{s}" for s in samples]
    done_fc = [x for x in done if x.startswith("epi2me_")]
    to_dos = [x for x in all_fc if x not in done_fc]

    dependencies = ":".join([f"${code}" for code in to_dos])

    if len(to_dos) > 0:
        print("To-Do: " + tool)
        with open(output + "/scripts/main.sh", "a") as f:
            f.write("\n# ont-methylDMR-kit")
            f.write(
                f"\nont_methyldmr_kit=$(sbatch --parsable --dependency=afterok:{dependencies} {job})\n"
            )
    else:
        if tool not in done:
            print("To-Do: " + tool)
            with open(output + "/scripts/main.sh", "a") as f:
                f.write("\n# ont-methylDMR-kit")
                f.write(f"\nont_methyldmr_kit=$(sbatch --parsable {job})\n")
        else:
            print("Done: " + tool)


def flair(toml_config, done):
    tool = "flair"

    output = toml_config["general"]["project_path"]
    email = toml_config["general"]["email"]
    genome = get_reference(toml_config["general"]["reference"])["fasta"]
    gtf = get_reference(toml_config["general"]["reference"])["gtf"]
    name = output.rstrip("/").split("/")[-2].split("_", 1)[1]
    username = os.environ.get("USER")

    for sample in toml_config["general"]["samples"]:
        # Make manifest (f"{sample}_manifest.txt"). Tab delimited file containing sample id, condition, batch,reads.fq, where reads.fq is the path to the sample fastq file.
        # TO-DO
        manifest = output + "/scripts/manifest_" + tool + "_" + sample + ".txt"
        with open(manifest, "w") as m:
            m.write(
                f"{sample}\tflair\t{name}\t/lustre10/scratch/{username}/{name}/results/flair/{sample}/{sample}.fastq"
            )

        bam = f"/lustre10/scratch/{username}/{name}/alignments/{sample}_sorted.bam"

        job = output + "/scripts/" + tool + "_" + sample + ".slurm"
        with open(
            TOOL_PATH
            + "main_pipelines/long-read/LongReadSequencingONT/template_flair.txt",
            "r",
        ) as f:
            slurm = f.read()
            slurm_filled = slurm.format(
                sample, email, name, output, genome, bam, gtf, manifest
            )

            with open(job, "w") as o:
                o.write(slurm_filled)

        samtools_name = f"samtools_{sample}"
        flair_name = f"flair_{sample}"

        if samtools_name not in done:
            print("To-Do: " + flair_name)
            with open(output + "/scripts/main.sh", "a") as f:
                f.write(f"\n# FLAIR for {sample}")
                f.write(
                    f"\n{flair_name}=$(sbatch --parsable --dependency=afterok:${samtools_name} {job})\n"
                )
        else:
            if flair_name not in done:
                print("To-Do: " + flair_name)
                with open(output + "/scripts/main.sh", "a") as f:
                    f.write(f"\n# FLAIR for {sample}")
                    f.write(f"\n{flair_name}=$(sbatch --parsable {job})\n")
            else:
                print("Done: " + flair_name)


def deepvariant(toml_config, done):
    tool = "deepvariant"

    output = toml_config["general"]["project_path"]
    email = toml_config["general"]["email"]
    genome = get_reference(toml_config["general"]["reference"])["fasta"]
    name = output.rstrip("/").split("/")[-2].split("_", 1)[1]
    username = os.environ.get("USER")

    for sample in toml_config["general"]["samples"]:
        bam = f"/lustre10/scratch/{username}/{name}/alignments/{sample}_sorted.bam"

        job = output + "/scripts/" + tool + "_" + sample + ".slurm"
        with open(
            TOOL_PATH
            + "main_pipelines/long-read/LongReadSequencingONT/template_deepvariant.txt",
            "r",
        ) as f:
            slurm = f.read()
            slurm_filled = slurm.format(sample, email, name, output, bam, genome)

            with open(job, "w") as o:
                o.write(slurm_filled)

        samtools_name = f"samtools_{sample}"
        deepvariant_name = f"deepvariant_{sample}"

        if samtools_name not in done:
            print("To-Do: " + deepvariant_name)
            with open(output + "/scripts/main.sh", "a") as f:
                f.write(f"\n# DeepVariant for {sample}")
                f.write(
                    f"\n{deepvariant_name}=$(sbatch --parsable --dependency=afterok:${samtools_name} {job})\n"
                )
        else:
            if deepvariant_name not in done:
                print("To-Do: " + deepvariant_name)
                with open(output + "/scripts/main.sh", "a") as f:
                    f.write(f"\n# DeepVariant for {sample}")
                    f.write(f"\n{deepvariant_name}=$(sbatch --parsable {job})\n")
            else:
                print("Done: " + deepvariant_name)


def cutesv(toml_config, done):
    tool = "cutesv"

    output = toml_config["general"]["project_path"]
    email = toml_config["general"]["email"]
    genome = get_reference(toml_config["general"]["reference"])["fasta"]
    name = output.rstrip("/").split("/")[-2].split("_", 1)[1]
    username = os.environ.get("USER")
    samples = toml_config["general"]["samples"]
    str_samples = " ".join(samples)

    job = output + "/scripts/" + tool + ".slurm"
    with open(
        TOOL_PATH
        + "main_pipelines/long-read/LongReadSequencingONT/template_cutesv.txt",
        "r",
    ) as f:
        slurm = f.read()
        slurm_filled = slurm.format(email, output, str_samples, username, name, genome)

        with open(job, "w") as o:
            o.write(slurm_filled)

    all_fc = [f"samtools_{s}" for s in samples]
    done_fc = [x for x in done if x.startswith("samtools")]
    to_dos = [x for x in all_fc if x not in done_fc]

    dependencies = ":".join([f"${code}" for code in to_dos])

    if len(to_dos) > 0:
        print("To-Do: " + tool)
        with open(output + "/scripts/main.sh", "a") as f:
            f.write("\n# cuteSV")
            f.write(
                f"\ncutesv=$(sbatch --parsable --dependency=afterok:{dependencies} {job})\n"
            )
    else:
        if tool not in done:
            print("To-Do: " + tool)
            with open(output + "/scripts/main.sh", "a") as f:
                f.write("\n# cuteSV")
                f.write(f"\ncutesv=$(sbatch --parsable {job})\n")
        else:
            print("Done: " + tool)


def transfer(toml_config, done):
    tool = "transfer"

    threads = "1"
    memory = "1"
    time = "00-01:00"

    output = toml_config["general"]["project_path"]
    email = toml_config["general"]["email"]
    name = output.rstrip("/").split("/")[-2].split("_", 1)[1]
    username = os.environ.get("USER")
    scratch = f"/lustre10/scratch/{username}/{name}"

    cmd = ["rsync", "-avxH", "--no-g", "--no-p", "--partial", scratch, output]

    command_str = "\n".join(cmd)

    # Make script but do not launch automatically
    create_script(tool, threads, memory, time, output, email, command_str, "")


def cleanup(toml_config, done):
    # Simple function to remove redundant files and cleanup structure
    tool = "cleanup"
    output = toml_config["general"]["project_path"]
    threads = "1"
    memory = "1"
    time = "00-01:00"
    email = toml_config["general"]["email"]
    flowcells = toml_config["general"]["fc_dir_names"]
    samples = toml_config["general"]["samples"]

    # Build cleanup commands
    commands = []

    # Move all logs to scripts/logs
    commands.append(f"mv {output}/*.log {output}/scripts/logs/")

    # Move all files in the main directory to scripts
    commands.append(
        f"mv -t {output}/scripts/ {output}/*.txt {output}/*.sh {output}/*.py {output}/*.slurm"
    )

    # Remove longreadsum useless output directory
    commands.append(f"rm -r {output}/output_LongReadSum")
    commands.append(f"rm -r {output}/log_output.log")

    # Remove empty dorado_demux.out
    commands.append(f"rm {output}/dorado_demux_run*.out")

    # Remove tmp from epi2me
    commands.append(f"rm -r {output}/work")
    commands.append(f"rm -r {output}/output")

    # Remove temp directories/files for each flowcell (MinKNOW-related OR redundant from merge)
    for flow in flowcells:
        commands.append(f"rm -r {output}/{flow}/fastq_*")
        commands.append(f"rm -r {output}/{flow}/bam_*")
        commands.append(f"rm -r {output}/{flow}/alignments/*.bam")
        commands.append(f"rm -r {output}/{flow}/alignments/*.bam.bai")
        commands.append(f"rm {output}/{flow}/main_reports/sequencing_summary*.txt")

    # Join all commands into a single string
    command_str = "\n".join(commands)

    # Create slurm job (Might be wasteful of server resources?)
    job = create_script(tool, threads, memory, time, output, email, command_str, "")

    # dependencies = ":".join([f"$epi2me_{sample}" for sample in samples])

    # with open(output + "/scripts/main.sh", "a") as f:
    #     f.write("\n# Cleanup temporary files and logs\n")
    #     f.write(f"\nsbatch --dependency=afterok:{dependencies} {job}\n")


# Launches main function
if __name__ == "__main__":
    main()
