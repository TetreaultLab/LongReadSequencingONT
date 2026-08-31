#!/bin/bash

CONFIG_FILE="$1"

# ------------------------------------------------------------------------------
# 1. Parse Project Name & Directory from TOML Config
# ------------------------------------------------------------------------------
if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: Configuration file '$CONFIG_FILE' not found." >&2
    exit 1
fi

PROJECT_PATH=$(grep -E '^\s*project_path\s*=' "$CONFIG_FILE" | cut -d'=' -f2 | tr -d ' "' | head -n 1)

ORIGINAL_SCRIPT="${PROJECT_PATH}/scripts/main.sh"
OUTPUT_SCRIPT="${PROJECT_PATH}/scripts/resubmit.sh"
DONE_FILE="${PROJECT_PATH}/scripts/steps_done.txt"
LOG_DIR="."

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

# Set log directory
if ls "${PROJECT_PATH}/scripts"/*.log >/dev/null 2>&1; then
    LOG_DIR="${PROJECT_PATH}/scripts"
else
    LOG_DIR="${PROJECT_PATH}"
fi

echo "Project Name: $PROJECT_NAME"

# ------------------------------------------------------------------------------
# 2. Initialize Output Script
# ------------------------------------------------------------------------------
cat << EOF > "$OUTPUT_SCRIPT"
#!/bin/bash
# Generated automatically on $(date)

CLEANUP_PATTERN="cleanup"
ACTIVE_CLEANUP_IDS=\$(squeue -u "\$USER" -h -o "%i %j" | awk -v pat="\$CLEANUP_PATTERN" '\$2 ~ pat {print \$1}')

if [ -n "\$ACTIVE_CLEANUP_IDS" ]; then
    for job_id in \$ACTIVE_CLEANUP_IDS; do
        echo "[CANCELING] Canceling active cleanup job: \$job_id"
        scancel "\$job_id"
    done
fi

EOF

# ------------------------------------------------------------------------------
# 3. Fetch Active Slurm Jobs (Using %200j to prevent name truncation!)
# ------------------------------------------------------------------------------
declare -A RUNNING_JOBS
while read -r job_id job_name; do
    if [ -n "$job_name" ]; then
        # Store exact job name
        RUNNING_JOBS["$job_name"]="$job_id"
    fi
done < <(squeue -u "$USER" -h -o "%i %200j")

# Helper functions
extract_core_name_from_log() {
    local log_file="$1"
    basename "$log_file" | sed -E 's/\.(rc|rg)[0-9]+\.[0-9]+\.log$//'
}

extract_job_id_from_log() {
    local log_file="$1"
    basename "$log_file" | awk -F'.' '{print $(NF-1)}'
}

# ------------------------------------------------------------------------------
# 4. Status Check Logic with Substring/Fuzzy Slurm Matching
# ------------------------------------------------------------------------------
check_status() {
    local tool_name="$1" # e.g., samtools_2C_CTL_SCD5

    # Check 1: steps_done.txt (Exact or substring match)
    if [ -f "$DONE_FILE" ] && grep -q -E "(^|[^a-zA-Z0-9_])${tool_name}([^a-zA-Z0-9_]|$)" "$DONE_FILE"; then
        echo "DONE"
        return
    fi

    # Check 2: Active in squeue (Direct lookup or partial name match)
    if [ -n "${RUNNING_JOBS[$tool_name]}" ]; then
        echo "RUNNING:${RUNNING_JOBS[$tool_name]}"
        return
    fi

    # Substring search in active jobs (handles custom SBATCH --job-name overrides)
    for active_name in "${!RUNNING_JOBS[@]}"; do
        if [[ "$active_name" == *"$tool_name"* ]] || [[ "$tool_name" == *"$active_name"* ]]; then
            echo "RUNNING:${RUNNING_JOBS[$active_name]}"
            return
        fi
    done

    # Check 3: Check log directory for matching logs
    local existing_log
    existing_log=$(ls -1t "${LOG_DIR}/${tool_name}"*.log 2>/dev/null | head -n 1)

    if [ -n "$existing_log" ]; then
        local log_core_name
        log_core_name=$(extract_core_name_from_log "$existing_log")

        if [[ "$log_core_name" == *"$tool_name"* ]]; then
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
# 5. Parse main.sh line by line
# ------------------------------------------------------------------------------
sanitize_line() {
    local line="$1"

    # Step 1: Replace multiple colons with a single colon
    line=$(echo "$line" | sed -E 's/:+/:/g')

    # Step 2: Remove a leading colon right after afterok: (e.g. afterok::job -> afterok:job)
    line=$(echo "$line" | sed -E 's/afterok::*/afterok:/g')

    # Step 3: Remove a trailing colon before a space (e.g. job: /path -> job /path)
    line=$(echo "$line" | sed -E 's/:[[:space:]]/ /g')

    # Step 4: If ALL dependencies were empty, remove the orphan flag completely
    line=$(echo "$line" | sed -E 's/--dependency=afterok:[[:space:]]/ /g')

    echo "$line"
}

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
                    echo "# [COMPLETED] $tool_name" >> "$OUTPUT_SCRIPT"
                    ;;
                "RUNNING")
                    running_id=$(echo "$status_info" | cut -d':' -f2)
                    echo "# [RUNNING] $tool_name (Job ID: $running_id)" >> "$OUTPUT_SCRIPT"
                    echo "DEPS+=(\"${running_id}\")" >> "$OUTPUT_SCRIPT"
                    ;;
                "NEED_RUN")
                    # Clean empty dependency slots before writing
                    cleaned_line=$(sanitize_line "$line")
                    echo "$cleaned_line" >> "$OUTPUT_SCRIPT"
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
                    echo "# [RUNNING] $tool_name (Job ID: $running_id)" >> "$OUTPUT_SCRIPT"
                    echo "${var_name}=\"${running_id}\"" >> "$OUTPUT_SCRIPT"
                    ;;
                "NEED_RUN")
                    # Clean empty dependency slots before writing
                    cleaned_line=$(sanitize_line "$line")
                    echo "$cleaned_line" >> "$OUTPUT_SCRIPT"
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
echo "Generated executable."
echo "Launch with 'bash $OUTPUT_SCRIPT'"
