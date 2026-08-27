#!/bin/bash

ORIGINAL_SCRIPT="scripts/main.sh"
OUTPUT_SCRIPT="scripts/main_resubmit.sh"
DONE_FILE="scripts/steps_done.txt"
LOG_DIR="."

# Initialize output script
echo "#!/bin/bash" > "$OUTPUT_SCRIPT"
echo "# Generated automatically on $(date)" >> "$OUTPUT_SCRIPT"
echo "" >> "$OUTPUT_SCRIPT"
echo "DEPS=()" >> "$OUTPUT_SCRIPT"
echo "" >> "$OUTPUT_SCRIPT"

# 1. Capture active Slurm jobs into an associative array: RUNNING_JOBS["job_name"]="job_id"
declare -A RUNNING_JOBS
while read -r job_id job_name; do
    if [ -n "$job_name" ]; then
        RUNNING_JOBS["$job_name"]="$job_id"
    fi
done < <(squeue -u "$USER" -h -o "%i %j")

# Helper: Extract core tool name from log filename pattern (stripping .rc/.rg host and jobid)
# Example: "dorado_basecaller_..._a881e64d.rg12902.17612157.log" -> "dorado_basecaller_..._a881e64d"
extract_core_name_from_log() {
    local log_path="$1"
    basename "$log_path" | sed -E 's/\.(rc|rg)[0-9]+\.[0-9]+\.log$//'
}

# Helper: Extract Job ID from log filename pattern
# Example: "dorado_basecaller_..._a881e64d.rg12902.17612157.log" -> "17612157"
extract_job_id_from_log() {
    local log_path="$1"
    basename "$log_path" | awk -F'.' '{print $(NF-1)}'
}

# Core status resolver function
check_status() {
    local tool_name="$1"

    # --- Check 1: steps_done.txt ---
    if [ -f "$DONE_FILE" ] && grep -q -F "$tool_name" "$DONE_FILE"; then
        echo "DONE"
        return
    fi

    # --- Check 2: Active in squeue ---
    for active_name in "${!RUNNING_JOBS[@]}"; do
        if [[ "$active_name" == *"$tool_name"* ]] || [[ "$tool_name" == *"$active_name"* ]]; then
            echo "RUNNING:${RUNNING_JOBS[$active_name]}"
            return
        fi
    done

    # --- Check 3: Log files in directory using regex matching ---
    local existing_log
    existing_log=$(ls -1 "${LOG_DIR}/${tool_name}"*.log 2>/dev/null | head -n 1)

    if [ -n "$existing_log" ]; then
        local log_core_name
        log_core_name=$(extract_core_name_from_log "$existing_log")

        # Verify exact base match
        if [ "$log_core_name" = "$tool_name" ]; then
            local job_id
            job_id=$(extract_job_id_from_log "$existing_log")

            if [[ "$job_id" =~ ^[0-9]+$ ]]; then
                # Check Slurm accounting history if job cleared squeue
                local state
                state=$(sacct -j "$job_id" -n -o State%15 2>/dev/null | head -n 1 | tr -d ' ')

                if [ "$state" = "COMPLETED" ]; then
                    echo "DONE"
                    return
                elif [[ "$state" == "RUNNING" || "$state" == "PENDING" ]]; then
                    echo "RUNNING:${job_id}"
                    return
                fi
            fi
        fi
    fi

    echo "NEED_RUN"
}

# 2. Process main.sh line by line
while IFS= read -r line || [ -n "$line" ]; do

    if [[ "$line" =~ sbatch.*\.slurm ]]; then

        # Extract .slurm file name and derive the tool basename
        script_file=$(echo "$line" | grep -oE '[^/]+\.slurm')
        tool_name="${script_file%.slurm}"

        status_info=$(check_status "$tool_name")
        status_type=$(echo "$status_info" | cut -d':' -f1)

        # Handle Array Addition: DEPS+=($(sbatch ...))
        if [[ "$line" =~ ^DEPS\+=\(\$\(sbatch ]]; then
            case "$status_type" in
                "DONE")
                    echo "# [COMPLETED] $tool_name" >> "$OUTPUT_SCRIPT"
                    ;;
                "RUNNING")
                    running_id=$(echo "$status_info" | cut -d':' -f2)
                    echo "# [ATTACHED] $tool_name (Running in Slurm - Job ID: $running_id)" >> "$OUTPUT_SCRIPT"
                    echo "DEPS+=(\"${running_id}\")" >> "$OUTPUT_SCRIPT"
                    ;;
                "NEED_RUN")
                    echo "$line" >> "$OUTPUT_SCRIPT"
                    ;;
            esac

        # Handle Variable Assignment: var_name=$(sbatch ...)
        elif [[ "$line" =~ ^([a-zA-Z0-9_]+)=\$\(sbatch ]]; then
            var_name="${BASH_REMATCH[1]}"
            case "$status_type" in
                "DONE")
                    echo "# [COMPLETED] $tool_name" >> "$OUTPUT_SCRIPT"
                    echo "${var_name}=\"\"" >> "$OUTPUT_SCRIPT"
                    ;;
                "RUNNING")
                    running_id=$(echo "$status_info" | cut -d':' -f2)
                    echo "# [ATTACHED] $tool_name (Running in Slurm - Job ID: $running_id)" >> "$OUTPUT_SCRIPT"
                    echo "${var_name}=\"${running_id}\"" >> "$OUTPUT_SCRIPT"
                    ;;
                "NEED_RUN")
                    echo "$line" >> "$OUTPUT_SCRIPT"
                    ;;
            esac
        else
            echo "$line" >> "$OUTPUT_SCRIPT"
        fi

    else
        # Pass through comments, empty lines, and control flow (if/else)
        echo "$line" >> "$OUTPUT_SCRIPT"
    fi

done < "$ORIGINAL_SCRIPT"

chmod +x "$OUTPUT_SCRIPT"
echo "Done! Generated executable: $OUTPUT_SCRIPT"
