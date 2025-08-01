#!/bin/bash
# Usage: bash run_statistics.sh <INPUT_DIRECTORY> [<TEST_NAME>] [<CORRECTION_METHOD>]

# Get the directory and optional test and correction arguments
directory="$1"
test_name="${2:-t-test}"           # Default to t-test if not provided
correction="${3:None}"             # Default to none if not provided

# Check if the directory exists
if [ ! -d "$directory" ]; then
    echo "Error: Directory not found."
    exit 1
fi

# Iterate over all pairs of files, excluding identical pairs
for file1 in "$directory"/*annotated.txt; do
    for file2 in "$directory"/*annotated.txt; do
        # Skip if the files are the same or if file1 is already compared with file2
        if [ "$file1" != "$file2" ]; then
            echo "Running statistics for: $file1 vs $file2 with $test_name and $correction"
            python scripts/fragpipe/statistics_multipletests.py \
                "$file1" "$file2" "$test_name" "$correction"
        fi
    done
done