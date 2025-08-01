#!/bin/bash
# Usage: bash run_volcano_list.sh <INPUT_DIRECTORY> 

# Get the directory and optional test and correction arguments
directory="$1"

# Check if the directory exists
if [ ! -d "$directory" ]; then
    echo "Error: Directory not found."
    exit 1
fi

# Iterate over all files
for file in "$directory"/*.txt; do
	echo "Creating volcano plot and significant list for "$file". "
        python ./scripts/fragpipe/volcano_plot.py "$file" figures/proteomics/
        python ./scripts/fragpipe/significant_list.py "$file" figures/proteomics/
done