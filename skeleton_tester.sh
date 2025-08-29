#!/bin/bash
set -euo pipefail

# Prompt for the source and destination directories
echo ">>> Please enter the source directory you want to copy for testing:"
read -r SRC

echo ">>> Please enter the destination directory - Important: Make sure to not use an existing directory!"
read -r DEST

# Prompt for dry-run option
echo ">>> Do you want to perform a dry-run (no actual copying) [y/n]?"
read -r DRY_ANSWER

if [[ "$DRY_ANSWER" =~ ^[Yy]$ ]]; then
    DRY_RUN=1
    echo ">>> Starting dry-run (no files will be copied)."
else
    DRY_RUN=0
    echo ">>> Starting actual skeleton copy."
fi

# Prompt for number of files
echo ">>> How many files per folder do you want to copy? (default = 3)"
read -r FILE_COUNT
if [[ -z "$FILE_COUNT" ]]; then
    FILE_COUNT=3
fi
echo ">>> Will copy up to $FILE_COUNT files per folder."

echo ">>> Source:      $SRC"
echo ">>> Destination: $DEST"
echo

# Normalize SRC to always have a trailing slash (so substring removal works)
SRC="${SRC%/}/"

# Make sure the destination directory doesn't exist already for a clean test run
if [[ -d "$DEST" ]]; then
    echo "Error: Destination directory $DEST already exists. Please provide a new, non-existing directory."
    exit 1
fi

mkdir -p "$DEST"

# Loop through each SamplePool
for samplepool in "$SRC"*/; do
    [[ -d "$samplepool" ]] || continue  # skip files
    rel_samplepool="${samplepool#$SRC}"

    if [[ $DRY_RUN -eq 1 ]]; then
        echo "    [DRY-RUN] Creating directory: $DEST/$rel_samplepool"
    else
        mkdir -p "$DEST/$rel_samplepool"
    fi

    # Every folder directly under SamplePool is a flowcell
    for flowcell in "$samplepool"*/; do
        [[ -d "$flowcell" ]] || continue
        rel_flowcell="${flowcell#$SRC}"
        dest_flowcell="$DEST/$rel_flowcell"

        echo ">>> Flowcell: $rel_flowcell"

        if [[ $DRY_RUN -eq 1 ]]; then
            echo "    [DRY-RUN] Creating directory: $dest_flowcell"
        else
            mkdir -p "$dest_flowcell"
        fi

        # Copy summary files at flowcell root (first level only)
        for file in "$flowcell"*; do
            if [[ -f "$file" ]] && [[ "$file" =~ \.(txt|csv|json|html|tsv|md)$ ]]; then
                relfile="${file#$SRC}"
                if [[ $DRY_RUN -eq 1 ]]; then
                    echo "    [DRY-RUN] Copying summary: $relfile"
                else
                    echo "    Copying summary: $relfile"
                    cp "$file" "$DEST/$relfile"
                fi
            fi
        done

        # Handle pod5 (search up to 2 levels deep inside flowcell)
        pod5_dirs=$(find "$flowcell" -maxdepth 2 -type d -name "pod5" 2>/dev/null || true)
        for pod5dir in $pod5_dirs; do
            rel_pod5="${pod5dir#$SRC}"
            if [[ $DRY_RUN -eq 1 ]]; then
                echo "    [DRY-RUN] Creating pod5 directory: $DEST/$rel_pod5"
            else
                mkdir -p "$DEST/$rel_pod5"
            fi
            files=$(ls -1 "$pod5dir" | sort | head -n "$FILE_COUNT" || true)
            for fname in $files; do
                file="$pod5dir/$fname"
                [[ -f "$file" ]] || continue
                relfile="${file#$SRC}"
                if [[ $DRY_RUN -eq 1 ]]; then
                    echo "    [DRY-RUN] Copying pod5 file: $relfile"
                else
                    echo "    Copying pod5 file: $relfile"
                    cp "$file" "$DEST/$relfile"
                fi
            done
        done

        # Handle fastq_pass and fastq_fail (files are in barcode subfolders)
        for subparent in fastq_pass fastq_fail; do
            parent="$flowcell/$subparent"
            [[ -d "$parent" ]] || continue
            for barcode in "$parent"/*; do
                [[ -d "$barcode" ]] || continue
                rel_barcode="${barcode#$SRC}"
                if [[ $DRY_RUN -eq 1 ]]; then
                    echo "    [DRY-RUN] Creating directory: $DEST/$rel_barcode"
                else
                    mkdir -p "$DEST/$rel_barcode"
                fi
                files=$(ls -1 "$barcode" | sort | head -n "$FILE_COUNT" || true)
                for fname in $files; do
                    file="$barcode/$fname"
                    [[ -f "$file" ]] || continue
                    relfile="${file#$SRC}"
                    if [[ $DRY_RUN -eq 1 ]]; then
                        echo "    [DRY-RUN] Copying $subparent file: $relfile"
                    else
                        echo "    Copying $subparent file: $relfile"
                        cp "$file" "$DEST/$relfile"
                    fi
                done
            done
        done

        # Handle other subdirectories fully
        for subdir in "$flowcell"*/; do
            [[ -d "$subdir" ]] || continue
            subname=$(basename "$subdir")
            [[ "$subname" =~ pod5|fastq_pass|fastq_fail|alignments ]] && continue
            rel_subdir="${subdir#$SRC}"

            if [[ $DRY_RUN -eq 1 ]]; then
                echo "    [DRY-RUN] Creating directory: $DEST/$rel_subdir"
            else
                mkdir -p "$DEST/$rel_subdir"
            fi

            for file in "$subdir"/*; do
                [[ -f "$file" ]] || continue
                relfile="${file#$SRC}"
                if [[ $DRY_RUN -eq 1 ]]; then
                    echo "    [DRY-RUN] Copying other subdir file: $relfile"
                else
                    echo "    Copying other subdir file: $relfile"
                    cp "$file" "$DEST/$relfile"
                fi
            done
        done
    done
done

echo
if [[ $DRY_RUN -eq 1 ]]; then
    echo ">>> Dry-run completed. No files were actually copied."
else
    echo ">>> Skeleton copy completed safely."
fi
