#!/bin/bash

# Create the output directory if it doesn't exist
if [ ! -d "3_fixed_pdb" ]; then
    mkdir "3_fixed_pdb"
fi

# Iterate over all aligned PDB files in the input directory
for in_file in ./2_aligned_pdb/*_aligned.pdb; do
    echo "Fixing $(basename "$in_file")"

    # Extract the basename and replace "_aligned" with "_fixed"
    output_file="3_fixed_pdb/$(basename "$in_file" _aligned.pdb)_fixed.pdb"

    # Run pdbfixer with the specified options
    pdbfixer "$in_file" \
        --output="$output_file" \
        --add-atoms=heavy \
        --keep-heterogens=all \
        --replace-nonstandard \
        --verbose
done