#!/bin/bash

# Save the original directory
orig_dir=$(pwd)

# Store the relative paths of PyMOL-APBS visualization scripts
isurf_pml_path="${orig_dir}/0_scripts/pymol_apbs_isosurface_vis.pml"

# Loop through all .pqr files in the specified directory
for pqr_file in ./4_pqr_apbs_files/*/*_fixed.pqr; do
    # Change to the directory containing the .pqr file and run PyMOL-APBS visualization scripts
    (
        cd "$(dirname "$pqr_file")" || exit
        pymol -c "$isurf_pml_path"
    )
done

# Return to the original directory
cd "$orig_dir"
