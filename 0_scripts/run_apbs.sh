#!/bin/bash

# Save the original directory
orig_dir=$(pwd)

# Loop through all .in files in the specified directory
for in_file in ./4_pqr_apbs_files/*/*.in; do
    # Change to the directory containing the .in file and run the apbs command
    (cd "$(dirname "$in_file")" && apbs "$(basename "$in_file")")
done

# Return to the original directory
cd "$orig_dir"