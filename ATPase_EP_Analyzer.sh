#!/bin/bash

echo "ATPase Electrostatic Potential Analyzer by Islam K. Matar"

# Activate Conda environment
eval "$(conda shell.bash hook)"
conda activate atpase_analyzer

# Step 1: Run rotor_ring_detection.py
python 0_scripts/rotor_ring_detection.py

# Step 2: Run PyMOL alignment pipeline
pymol -c ATPase_Alignment_Pipeline.pml

# Step 3: Remove unnecessary file
rm UNK.cif

# Step 4: Run pdbfixer on aligned ATPases PDB files
bash 0_scripts/run_pdbfixer.sh

# Step 5: Run custom pdb2pqr python script
python 0_scripts/pdb2pqr.py

# Step 6: Run APBS calculations
export OMP_NUM_THREADS=50
bash 0_scripts/run_apbs.sh

# Step 7: Perform CT-scan downstream analysis
python 0_scripts/apbs_ctscan_analysis.py

# Step 8: Generate isosurfaces using PyMOL
bash 0_scripts/pymol_apbs_analysis.sh

# Step 9: Compress APBS output files for download
tar -czvf 4_pqr_apbs_files.tar.gz 4_pqr_apbs_files/

# Step 10: Compute file size before download
du -sh 4_pqr_apbs_files.tar.gz

echo "Pipeline execution completed."
