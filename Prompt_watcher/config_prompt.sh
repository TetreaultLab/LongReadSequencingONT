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

# Project path
echo
echo ">>> Project folder name (This includes the name you gave the sample; i.e. 2025-08-25_FXN-WGS-Batch1/FXN_N12)"
project_path=$(prompt "project_path" "Project path" "")
echo "    project_path = \"/lustre09/project/6019267/shared/projects/Nanopore_Dock/$project_path\"" >> "$CONFIG"

# Barcode
barcode=$(prompt "barcode" "Barcode (first number or leave empty for none)" "")
echo "    barcode = $barcode" >> "$CONFIG"

# Samples
echo
echo ">>> Samples: Enter sample names separated by space. Ex: Control1 Case1"
read -p "> " samples
samples_list="[\"$(echo $samples | sed 's/ /\",\"/g')\"]"
echo "    samples = $samples_list" >> "$CONFIG"

# Conditions
echo
echo ">>> Conditions: Enter conditions matching the samples, space-separated. Ex: Ctrl Disease"
read -p "> " conditions
conditions_list="[\"$(echo $conditions | sed 's/ /\",\"/g')\"]"
echo "    conditions = $conditions_list" >> "$CONFIG"

# Seq type
seq_type=$(prompt "seq_type" "Seq type (WGS, RNA, Targeted)" "")
echo "    seq_type = \"$seq_type\"" >> "$CONFIG"

# Kit
echo
echo ">>> Kit options:"
echo "    WGS: [\"SQK-RBK114-24\", \"SQK-NBD114-24\", \"SQK-LSK114\"]"
echo "    RNA: [\"SQK-PCB114-24\"]"
echo "    Targeted: [\"SQK-NBD114-24\"]"
kit=$(prompt "kit" "Kit (exact option)" "")
echo "    kit = \"$kit\"" >> "$CONFIG"

# Email
echo
echo ">>> Your e-mail address to receive notification"
email=$(prompt "email" "Email" "")
echo "    email = \"$email\"" >> "$CONFIG"

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
