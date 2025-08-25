#!/bin/bash

WATCH_DIRS=("/data" "/seqdata")
SCRIPT="/home/grub/Desktop/config_prompt.sh"

declare -A WATCHED_PROJECTS

watch_project() {
    local PROJECT="$1"
    if [[ -z "${WATCHED_PROJECTS[$PROJECT]}" ]]; then
        WATCHED_PROJECTS["$PROJECT"]=1
        echo "Watching project folder: $PROJECT"

        # Watch this project folder for new directories
        inotifywait -m -e create -e moved_to --format '%w%f' "$PROJECT" | while read NEWPATH
        do
            # Make sure it is a directory (could get a file event)
            if [[ -d "$NEWPATH" ]]; then
                echo "New subfolder detected: $NEWPATH"
                if [ ! -f "$NEWPATH/.watcher_seen" ]; then
                    touch "$NEWPATH/.watcher_seen"
                    gnome-terminal --working-directory="$NEWPATH" -- bash -c "$SCRIPT '$NEWPATH'; exec bash"
                fi
            fi
        done &
    fi
}

# Watch top-level directories for new project folders
inotifywait -m -e create -e moved_to --format '%w%f' "${WATCH_DIRS[@]}" | while read NEWPATH
do
    BASENAME=$(basename "$NEWPATH")
    if [[ -d "$NEWPATH" && "$BASENAME" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}.* ]]; then
        echo "New project folder detected: $NEWPATH"

        # 1) Prompt for any subfolders created at the same time as the project
        for SUB in "$NEWPATH"/*/; do
            [[ -d "$SUB" ]] || continue
            if [ ! -f "$SUB/.watcher_seen" ]; then
                touch "$SUB/.watcher_seen"
                echo "Subfolder found inside new project: $SUB"
                gnome-terminal --working-directory="$SUB" -- bash -c "$SCRIPT '$SUB'; exec bash"
            fi
        done

        # 2) Start watching this new project
        watch_project "$NEWPATH"
    fi
done &

# At startup, watch all existing project folders
for DIR in "${WATCH_DIRS[@]}"; do
    for PROJECT in "$DIR"/*; do
        if [[ -d "$PROJECT" && $(basename "$PROJECT") =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}.* ]]; then
            watch_project "$PROJECT"
        fi
    done
done

wait
