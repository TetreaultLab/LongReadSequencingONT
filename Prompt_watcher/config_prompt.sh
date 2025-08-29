#!/bin/bash
# Interactive config generator for Nanopore pipeline

TARGET_DIR="${1:-$(pwd)}"
CONFIG="$TARGET_DIR/config_initial.txt"
echo "[general]" > "$CONFIG"

prompt() {
	local varname=$1
    	local comment=$2
    	local default=$3
    	if [[ -n "$default" ]]; then
        	read -p "$comment (press Enter to use default: $default): " value
        	value=${value:-$default}
    	else
        	read -p "$comment: " value
    	fi
    	echo "$value"
}

# Pipeline valid options - To change if new options are added / options are removed : 
valid_kits=("SQK-RBK114-24" "SQK-NBD114-24" "SQK-LSK114" "SQK-PCB114-24") # Kits used in the lab
valid_seq=("WGS" "RNA" "Targeted") # Sequencing types handled by pipeline

# Project path
while true; do
	echo
    	echo ">>> Project folder name (This includes the name you gave the sample; i.e. 2025-08-25_FXN-WGS-Batch1/FXN_N12)"
    	project_path=$(prompt "project_path" "Project path" "")

    	# Must contain exactly one slash that separates project name from sample name
    	if [[ $(grep -o "/" <<< "$project_path" | wc -l) -ne 1 ]]; then
        	echo "⚠️  Warning: The path you've entered is in the wrong format. It should be 'ProjectName/SampleName'. Please try again."
        	continue
    	fi

    	# Extract the two informations 
    	project_prefix="${project_path%%/*}"  # Everything before the slash
    	project_suffix="${project_path##*/}"  # Everything after the slash

    	# Validate first part matches YYYY-MM-DD_XXX
    	if [[ ! "$project_prefix" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}_.{3,}$ ]]; then
        	echo "⚠️  Warning: Your project name \"$project_prefix\" is not correctly set for the pipeline to work."
        	echo "    Expected format: YYYY-MM-DD_ProjectName (e.g., 2025-08-25_FXN-WGS-Batch1)"
        	echo "    Please retry. If error persists, exit with CTRL + C."
		echo
		echo "    Important: You CANNOT rename the folder while MinKNOW is actively sequencing."
		echo "    Either" 
		echo "           1) Stop the sequencing + Restart with correct name + transfer pod5s"
		echo "           2) Rename after the run (may break reports)."
        	continue
    	fi

    	# When both checks are okay, we move on
    	echo "    project_path = \"/lustre09/project/6019267/shared/projects/Nanopore_Dock/$project_path\"" >> "$CONFIG"
    	break
done

# Samples
echo
echo ">>> Samples: Enter sample names separated by space. Ex: Control1 Case1"
read -p "> " samples
samples_array=($samples) 		# Make a an array to count samples
sample_count=${#samples_array[@]}	# Store count to compare with barcodes/conditions
samples_list="[\"$(echo $samples | sed 's/ /\",\"/g')\"]"
echo "    samples = $samples_list" >> "$CONFIG"

# Conditions
echo
echo ">>> Conditions: Enter conditions matching the samples, space-separated. Ex: Ctrl Disease"
while true; do
	read -p "> " conditions
    	conditions_array=($conditions) 		# Make an array just like samples
    	condition_count=${#conditions_array[@]}	# Store length to compare to samples

    	# Comparison loop
    	if [[ $condition_count -ne $sample_count ]]; then
        	echo "⚠️  Number of conditions ($condition_count) doesn't match the number of samples ($sample_count). Please try again."
        	continue
    	fi

    	# Format to list
    	conditions_list="[\"$(echo $conditions | sed 's/ /\",\"/g')\"]"
    	echo "    conditions = $conditions_list" >> "$CONFIG"
    	break
done

# Barcode
echo ">>> If you have consecutive barcodes, just input the number of the first one"
echo ">>> For non-consecutive barcodes, you need to input a space-separated list of all barcodes used in the order of samples"
while true; do
	echo
	barcode_input=$(prompt "barcode" "Barcode (Either single i.e. '1' or multiple i.e. '2, 5, 7, 3'; leave empty for none)" "")

	# If the input from the user is empty (no barcodes)
	if [[ -z "$barcode_input" ]]; then
        	echo "    barcode = []" >> "$CONFIG"
        	break
    	fi

    	# Split input into array
    	IFS=' ' read -r -a barcode_array <<< "$barcode_input"
	if [[ ${#barcode_array[@]} -eq 1 ]]; then
        	# Single barcode input — treated as starting barcode - no list format
        	echo "    barcode = ${barcode_array[0]}" >> "$CONFIG"
        	break
    	else
        	# Multiple barcode inputs — must match sample count - list format
        	if [[ ${#barcode_array[@]} -ne $sample_count ]]; then
            		echo "⚠️  Number of barcodes (${#barcode_array[@]}) doesn't match the number of samples ($sample_count). Please try again."
            		continue
        	fi

        	# Reformat array in list format for launch_pipeline.py
        	barcode_list="[${barcode_array[*]}]"
        	barcode_list=$(echo "$barcode_list" | sed 's/ /, /g')
        	echo "    barcode = $barcode_list" >> "$CONFIG"
        	break
    	fi
done

# Seq type
echo
echo ">>> Seq options:"
echo "    WGS (Whole-genome sequencing)"
echo "    RNA (Transcriptomic sequencing)"
echo "    Targeted (PCR amplified or other)"
while true; do
	echo
	seq_type=$(prompt "seq_type" "Seq type (WGS, RNA, Targeted)" "")
	
	# Check if input is in valid_kits array
        if [[ " ${valid_seq[@]} " =~ " $seq_type " ]]; then
		echo "    seq_type = \"$seq_type\"" >> "$CONFIG"
		break
	else
		echo "⚠️ Warning: The type you've entered is invalid. Please re-try."
	fi
done

# Kit
echo
echo ">>> Kit options:"
echo "    WGS: [\"SQK-RBK114-24\", \"SQK-NBD114-24\", \"SQK-LSK114\"]"
echo "    RNA: [\"SQK-PCB114-24\"]"
echo "    Targeted: [\"SQK-NBD114-24\"]"
while true; do
    	echo
	kit=$(prompt "kit" "Kit (exact option)" "")

    	# Check if input is in valid_kits array
    	if [[ " ${valid_kits[@]} " =~ " $kit " ]]; then
        	echo "    kit = \"$kit\"" >> "$CONFIG"
        	break
    	else
        	echo "⚠️ Warning: The kit name you've entered is invalid. Please re-try."
    	fi
done

# Email
echo 
echo ">>> Your e-mail address to receive notification"
while true; do
	echo
    	email=$(prompt "email" "Email" "")

    	# Simple check that the address somewhat looks valid - doesn't allow a skip
    	if [[ "$email" =~ ^[^@]+@[^@]+\.[^@]+$ ]]; then
        	echo "    email = \"$email\"" >> "$CONFIG"
        	break
    	else
        	echo "⚠️ Sorry but this doesn't seem like an email address. Please enter a valid email address (i.e. user@example.com)."
    	fi
done

# Analysis
echo
echo ">>> Analysis: Enter space-separated analysis modules. Leave empty to skip (default: [])."
echo "    If default is used, will only run Dorado & QC for basic alignment and output."
echo "    Options: none, splicing, methylation, SNP, SV, polya, repeats"
read -p "> " analysis
analysis_list="[\"$(echo $analysis | sed 's/ /\",\"/g')\"]"
echo "    analysis = $analysis_list" >> "$CONFIG"

# File type
echo
echo ">>> File type: pod5 is the only available option for now."
file_type=$(prompt "file_type" "File type" "pod5")
echo "    file_type = \"$file_type\"" >> "$CONFIG"

# Reference
echo
echo ">>> Reference genome: grch38 is the only available option for now."
reference=$(prompt "reference" "Reference" "grch38")
echo "    reference = \"$reference\"" >> "$CONFIG"

echo
echo ">>> Config saved to $CONFIG"
echo ">>> Now you need to use the script that creates all config files"
echo
echo "1. Activate the environment: source /lustre09/project/6019267/shared/tools/main_pipelines/long-read/launch_pipeline_env/bin/activate"
echo "(This steps activates all the pre-requisite modules for the python scripts)"
echo
echo "2. Launch the script: python3 /lustre09/project/6019267/shared/tools/main_pipelines/long-read/LongReadSequencingONT/launch_pipeline.py config_initial.txt" 
echo "This should automatically create and launch the pipeline :)" 
echo "If it doesn't, the pipeline can be launched manually: sh ./scripts/main.sh"
echo
echo
