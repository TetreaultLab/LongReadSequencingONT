#!/bin/bash

ORIGINAL_SCRIPT="scripts/main.sh"
OUTPUT_SCRIPT="scripts/resubmit.sh"
DONE_FILE="scripts/steps_done.txt"
LOG_DIR="."
CONFIG_FILE="config_final.toml"

# ------------------------------------------------------------------------------
# 1. Parse Project Name & Directory from TOML Config
# ------------------------------------------------------------------------------
if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: Configuration file '$CONFIG_FILE' not found." >&2
    exit 1
fi

# Extract raw project_path line
raw_path=$(grep -E '^\s*project_path\s*=' "$CONFIG_FILE" | sed -E 's/.*"([^"]+)".*/\1/')

if [ -z "$raw_path" ]; then
    echo "Error: Could not extract project_path from $CONFIG_FILE" >&2
    exit 1
fi

# Process path: remove trailing slash -> extract parent dir -> split date -> replace - with _
clean_path="${raw_path%/}"
parent_dir="${clean_path%/*}"
dir_name="${parent_dir##*/}"
PROJECT_NAME="${dir_name#*_}"

# Set log directory (uses the target script path)
LOG_DIR="${clean_path}/scripts"
[ ! -d "$LOG_DIR" ] && LOG_DIR="."

echo "[CONFIG] Project Name: $PROJECT_NAME"
echo "[CONFIG] Log Directory: $LOG_DIR"

# ------------------------------------------------------------------------------
# 2. Initialize output script & insert Automatic Cleanup Cancellation
# ------------------------------------------------------------------------------
cat << EOF > "$OUTPUT_SCRIPT"
#!/bin/bash
# Generated automatically on $(date)

# ------------------------------------------------------------------------------
# Cancel existing cleanup job if running/queued to prevent premature file deletion
# ------------------------------------------------------------------------------
CLEANUP_PATTERN="cleanup"
ACTIVE_CLEANUP_IDS=\$(squeue -u "\$USER" -h -o "%i %j" | awk -v pat="\$CLEANUP_PATTERN" '\$2 ~ pat {print \$1}')

if [ -n "\$ACTIVE_CLEANUP_IDS" ]; then
    for job_id in \$ACTIVE_CLEANUP_IDS; do
        echo "[CANCELING] Canceling active cleanup job: \$job_id"
        scancel "\$job_id"
    done
fi

DEPS=()

EOF

# ------------------------------------------------------------------------------
# 3. Fetch active Slurm jobs for current user
# ------------------------------------------------------------------------------
declare -A RUNNING_JOBS
while read -r job_id job_name; do
    [ -n "$job_name" ] && RUNNING_JOBS["$job_name"]="$job_id"
done < <(squeue -u "$USER" -h -o "%i %j")

# Helper functions for log processing
extract_core_name_from_log() {
    local log_file="$1"
    basename "$log_file" | sed -E 's/\.(rc|rg)[0-9]+\.[0-9]+\.log$//'
}

extract_job_id_from_log() {
    local log_file="$1"
    basename "$log_file" | awk -F'.' '{print $(NF-1)}'
}

# ------------------------------------------------------------------------------
# 4. Status Check Logic
# ------------------------------------------------------------------------------
check_status() {
    local tool_name="$1"

    # Check 1: steps_done.txt (Exact line match)
    if [ -f "$DONE_FILE" ] && grep -q -x "$tool_name" "$DONE_FILE"; then
        echo "DONE"
        return
    fi

    # Check 2: Active in squeue (Exact Job Name match)
    if [ -n "${RUNNING_JOBS[$tool_name]}" ]; then
        echo "RUNNING:${RUNNING_JOBS[$tool_name]}"
        return
    fi

    # Check 3: Check logs in scoped LOG_DIR
    local existing_log
    existing_log=$(ls -1t "${LOG_DIR}/${tool_name}".*.log 2>/dev/null | head -n 1)

    if [ -n "$existing_log" ]; then
        local log_core_name
        log_core_name=$(extract_core_name_from_log "$existing_log")

        if [ "$log_core_name" = "$tool_name" ]; then
            local job_id
            job_id=$(extract_job_id_from_log "$existing_log")

            if [[ "$job_id" =~ ^[0-9]+$ ]]; then
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

# ------------------------------------------------------------------------------
# 5. Parse main.sh and generate main_resubmit.sh
# ------------------------------------------------------------------------------
while IFS= read -r line || [ -n "$line" ]; do

    if [[ "$line" =~ sbatch.*\.slurm ]]; then

        script_file=$(echo "$line" | grep -oE '[^/]+\.slurm')
        tool_name="${script_file%.slurm}"

        status_info=$(check_status "$tool_name")
        status_type=$(echo "$status_info" | cut -d':' -f1)

        # Handle Array Addition: DEPS+=($(sbatch ...))
        if [[ "$line" =~ ^DEPS\+=\(\$\(sbatch ]]; then
            case "$status_type" in
                "DONE")
                    echo "# [SKIPPED] $tool_name (Completed)" >> "$OUTPUT_SCRIPT"
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
                    echo "# [SKIPPED] $tool_name (Completed)" >> "$OUTPUT_SCRIPT"
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
        echo "$line" >> "$OUTPUT_SCRIPT"
    fi

done < "$ORIGINAL_SCRIPT"

chmod +x "$OUTPUT_SCRIPT"
echo "Done! Generated executable. Launch with 'bash $OUTPUT_SCRIPT'"
