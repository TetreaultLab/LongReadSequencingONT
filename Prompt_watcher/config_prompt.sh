#!/bin/bash

# Interactive config generator for Nanopore pipeline

TARGET_DIR="${1:-$(pwd)}"
CONFIG="$TARGET_DIR/config_initial.txt"
echo "[general]" > "$CONFIG"

# Definition of inputs format for file configuration
prompt() {
	local varname=$1
	local comment=$2
	local default=$3
	if [[ -n "$default" ]]; then
		read -e -p "$comment (press Enter to use default: $default): " value
		value=${value:-$default}	
	else
		read -e -p "$comment: " value
	fi
	echo "$value"
}

# Pipeline valid options - Will need to add later
valid_kits=("SQK-RBK114-24" "SQK-NBD114-24" "SQK-LSK114" "SQK-PCB114-24")
valid_seq=("WGS" "RNA" "Targeted")

############################################
# Additional navigation added for base users
############################################
# Step counter to go back between them
current_step=1
total_steps=10

echo
echo "Base users guide:"
echo "  - Added option to use arrow keys to edit current line"
echo "  - Can now type 'back' to return to previous section"
echo

while [[ $current_step -le $total_steps ]]; do

case $current_step in

############################################
# 1) PROJECT PATH
############################################
1)
while true; do
	echo
	echo ">>> Project folder name (This includes the name you gave the sample; i.e. 2025-08-25_FXN-WGS-Batch1/FXN_N12)"
	read -e -p "Project path: " project_path

	if [[ "$project_path" == "back" ]]; then
		((current_step--))
		break
	fi

	# Format check
	if [[ $(grep -o "/" <<< "$project_path" | wc -l) -ne 1 ]]; then
		echo "⚠️  Warning: The path you've entered is in the wrong format. It should be 'ProjectName/SampleName'. Please try again."
		continue
	fi

	project_prefix="${project_path%%/*}"
	project_suffix="${project_path##*/}"

	if [[ ! "$project_prefix" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}_.{3,}$ ]]; then
		echo "⚠️ Warning: Your project name \"$project_prefix\" is not correctly set for the pipeline to work."
		echo "    Expected format: YYYY-MM-DD_ProjectName (e.g., 2025-08-25_FXN-WGS-Batch1)"
		echo "    Please retry. If error persists, exit with CTRL + C."
		echo
		echo "    Important: You CANNOT rename the folder while MinKNOW is actively sequencing."
		echo "    Either"
		echo "	    1) Stop the sequencing + Restart with correct name + transfer pod5s"
		echo "	    2) Rename after the run (may break reports)."
		continue
	fi

	echo "    project_path = \"/lustre09/project/6019267/shared/projects/Nanopore_Dock/$project_path\"" >> "$CONFIG"
	((current_step++))
	break
done
;;

############################################
# 2) SAMPLE NAMES
############################################
2)
echo
echo ">>> Samples: Enter sample names separated by space. Ex: Control1 Case1 UTMAB-43"
read -e -p "> " samples

if [[ "$samples" == "back" ]]; then
	((current_step--))
	continue
fi

samples_array=($samples)
sample_count=${#samples_array[@]}

# --- Duplicate warning (non-blocking) ---
declare -A seen
for s in "${samples_array[@]}"; do
	if [[ -n "${seen[$s]}" ]]; then
		echo "⚠️ WARNING: Duplicate sample detected -> '$s', go back if this is unexpected, proceed if desired."
	else
		seen[$s]=1
	fi
done
# --- Not adding further filters yet -----

samples_list="[\"$(echo $samples | sed 's/ /\",\"/g')\"]"
echo "    samples = $samples_list" >> "$CONFIG"
((current_step++))
;;

############################################
# 3) CONDITIONS
############################################
3)
echo
echo ">>> Conditions: Enter conditions matching the samples, space-separated. Ex: Ctrl Disease Ctrl"
echo ">>> For this metadata, there can be only a two/few groups, or all conditions can be different."
while true; do
	read -e -p "> " conditions

	if [[ "$conditions" == "back" ]]; then
		((current_step--))
		break
	fi

	conditions_array=($conditions)
	condition_count=${#conditions_array[@]}

	# Need to have as many conditions as there are samples
	if [[ $condition_count -ne $sample_count ]]; then
		echo "⚠️  Number of conditions ($condition_count) doesn't match the number of samples ($sample_count). Please try again."
		continue
	fi

	conditions_list="[\"$(echo $conditions | sed 's/ /\",\"/g')\"]"
	echo "    conditions = $conditions_list" >> "$CONFIG"
	((current_step++))
	break
done
;;

############################################
# 4) BARCODE
############################################
4)
echo ">>> If you have consecutive barcodes, just input the number of the first one (i.e. NB03 to NB011 --> '3')"
echo ">>> For non-consecutive barcodes, you need to input a space-separated list of all barcodes used in the order of samples"
echo ">>> The list should then look like : '2 5 7 3'. If no barcodes are used (i.e. LSK), leave empty."
while true; do
	echo
	read -e -p "Barcode: " barcode_input

	if [[ "$barcode_input" == "back" ]]; then
		((current_step--))
		break
	fi

	# Format check (otherwise allows user to input uninformative blocks)
	if [[ ! "$barcode_input" =~ ^[0-9[:space:]]*$ ]]; then
	echo "⚠️ Warning: Barcode input contains invalid characters. Only numbers and spaces are allowed."
	continue
	fi

	if [[ -z "$barcode_input" ]]; then
		echo "    barcode = []" >> "$CONFIG"
		((current_step++))
		break
	fi

	IFS=' ' read -r -a barcode_array <<< "$barcode_input"

	# Allow for cases 1) No barcode or 2) Only 1 barcode for > 2 samples
	if [[ ${#barcode_array[@]} -eq 1 ]]; then
		echo "    barcode = ${barcode_array[0]}" >> "$CONFIG"
		((current_step++))
		break
	fi

	# For non-consecutive barcodes, check that quantity matches number of samples
	if [[ ${#barcode_array[@]} -ne $sample_count ]]; then
		echo "⚠️  Number of barcodes (${#barcode_array[@]}) doesn't match the number of samples ($sample_count). Please try again."
		continue
	fi

	barcode_list="[${barcode_array[*]}]"
	barcode_list=$(echo "$barcode_list" | sed 's/ /, /g')
	echo "    barcode = $barcode_list" >> "$CONFIG"
	((current_step++))
	break
done
;;

############################################
# 5) SEQUENCING TYPE
############################################
5)
echo
echo ">>> Seq options:"
echo "    WGS (Whole-genome sequencing)"
echo "    RNA (Transcriptomic sequencing)"
echo "    Targeted (PCR amplified or other)"
while true; do
	echo
	read -e -p "Seq type: " seq_type

	if [[ "$seq_type" == "back" ]]; then
		((current_step--))
		break
	fi

	# Validation from modifyable list
	if [[ " ${valid_seq[@]} " =~ " $seq_type " ]]; then
		echo "    seq_type = \"$seq_type\"" >> "$CONFIG"
		((current_step++))
		break
	else
		echo "⚠️ Warning: The type you've entered is invalid. Please re-try."
	fi
done
;;

############################################
# 6) KIT NAME
############################################
6)
echo
echo ">>> Kit options:"
echo "    WGS: [\"SQK-RBK114-24\", \"SQK-NBD114-24\", \"SQK-LSK114\"]"
echo "    RNA: [\"SQK-PCB114-24\"]"
echo "    Targeted: [\"SQK-NBD114-24\"]"
while true; do
	echo
	read -e -p "Kit: " kit

	if [[ "$kit" == "back" ]]; then
		((current_step--))
		break
	fi

	# Consider adding a kit match to sequencing type? (May be redundant in pipeline)

	# Validation from modifyable list
	if [[ " ${valid_kits[@]} " =~ " $kit " ]]; then
		echo "    kit = \"$kit\"" >> "$CONFIG"
		((current_step++))
		break
	else
		echo "⚠️ Warning: The kit name you've entered is invalid. Please re-try."
	fi
done
;;

############################################
# 7) EMAIL ADDRESS
############################################
7)
echo
echo ">>> Your e-mail address to receive notification"
while true; do
	echo
	read -e -p "Email: " email

	if [[ "$email" == "back" ]]; then
		((current_step--))
		break
	fi

	# Quick check if looks like an email address
	if [[ "$email" =~ ^[^@]+@[^@]+\.[^@]+$ ]]; then
		echo "    email = \"$email\"" >> "$CONFIG"
		((current_step++))
		break
	else
		echo "⚠️ Sorry but this doesn't seem like an email address. Please enter a valid email address (i.e. user@example.com)."
	fi
done
;;

############################################
# 8) ANALYSIS
############################################
8)
echo
echo ">>> Analysis: Enter space-separated analysis modules. Leave empty to skip (default: [])."
echo "    If default is used, will only run Dorado & QC for basic alignment and output."
echo "    Options: splicing, methylation, SNP, SV, CNV, polya, repeats, phasing"
read -e -p "> " analysis

if [[ "$analysis" == "back" ]]; then
	((current_step--))
	continue
fi

analysis_list="[\"$(echo $analysis | sed 's/ /\",\"/g')\"]"
echo "    analysis = $analysis_list" >> "$CONFIG"
((current_step++))
;;

############################################
# 9) FILE TYPE
############################################
9)
echo
echo ">>> File type: pod5 is the only available option for now."
read -e -p "File type (press Enter to use default: pod5): " file_type
file_type=${file_type:-pod5}

if [[ "$file_type" == "back" ]]; then
	((current_step--))
	continue
fi

echo "    file_type = \"$file_type\"" >> "$CONFIG"
((current_step++))
;;

############################################
# 10) REFERENCE
############################################
10)
echo
echo ">>> Reference genome: grch38 is the only available option for now."
read -e -p "Reference (press Enter to use default: grch38): " reference
reference=${reference:-grch38}

if [[ "$reference" == "back" ]]; then
	((current_step--))
	continue
fi

echo "    reference = \"$reference\"" >> "$CONFIG"
((current_step++))
;;

esac
done

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

############################################
# FINAL REVIEW / RE-EDIT MENU
############################################

while true; do
	echo
	echo "=========== FINAL REVIEW ==========="
	echo "1) Proj. path	: $project_path"
	echo "2) Samples	: $samples"
	echo "3) Condition	: $conditions"
	echo "4) Barcode	: $barcode_input"
	echo "5) Seq type	: $seq_type"
	echo "6) Kit	: $kit"
	echo "7) Email	: $email"
	echo "8) Analysis	: $analysis"
	echo "9) File type	: $file_type"
	echo "10) Reference	: $reference"
	echo "------------------------------------"
	echo "0) Finish and keep config"
	echo

	read -e -p "Select number to re-edit (0 to finish): " choice

	# Jump back into a specific section (not throughout the whole thing)
	case $choice in
	1)
		echo
		read -e -p "Re-enter Project path: " project_path
		sed -i "s|project_path = .*|project_path = \"/lustre09/project/6019267/shared/projects/Nanopore_Dock/$project_path\"|" "$CONFIG"
		;;
	2)
		echo
		read -e -p "Re-enter Samples: " samples
		samples_list="[\"$(echo $samples | sed 's/ /\",\"/g')\"]"
		sed -i "s|samples = .*|samples = $samples_list|" "$CONFIG"
		;;
	3)
		echo
		read -e -p "Re-enter Conditions: " conditions
		conditions_list="[\"$(echo $conditions | sed 's/ /\",\"/g')\"]"
		sed -i "s|conditions = .*|conditions = $conditions_list|" "$CONFIG"
		;;
	4)
		echo
		read -e -p "Re-enter Barcode: " barcode_input
		if [[ -z "$barcode_input" ]]; then
			barcode_value="[]"
		else
		
			IFS=' ' read -r -a barcode_array <<< "$barcode_input"
			if [[ ${#barcode_array[@]} -eq 1 ]]; then
				barcode_value="${barcode_array[0]}"
			else
				barcode_value="[${barcode_array[*]}]"
				barcode_value=$(echo "$barcode_value" | sed 's/ /, /g')
			fi
		fi
		sed -i "s|barcode = .*|barcode = $barcode_value|" "$CONFIG"
		;;
	5)
		read -e -p "Re-enter Seq type: " seq_type
		sed -i "s|seq_type = .*|seq_type = \"$seq_type\"|" "$CONFIG"
		;;
	6)
		read -e -p "Re-enter Kit: " kit
		sed -i "s|kit = .*|kit = \"$kit\"|" "$CONFIG"
		;;
	7)
		read -e -p "Re-enter Email: " email
		sed -i "s|email = .*|email = \"$email\"|" "$CONFIG"
		;;
	8)
		read -e -p "Re-enter Analysis: " analysis
		analysis_list="[\"$(echo $analysis | sed 's/ /\",\"/g')\"]"
		sed -i "s|analysis = .*|analysis = $analysis_list|" "$CONFIG"
		;;
	9)
		read -e -p "Re-enter File type: " file_type
		sed -i "s|file_type = .*|file_type = \"$file_type\"|" "$CONFIG"
		;;
	10)
		read -e -p "Re-enter Reference: " reference
		sed -i "s|reference = .*|reference = \"$reference\"|" "$CONFIG"
		;;
	0)
		echo
		echo "Config finalized."
		break
		;;
	# If somehow anything else
	*)
		echo "Invalid choice."
		;;
	esac
done
