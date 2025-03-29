# Loading dependencies and custom functions
run ./0_scripts/centroid.py
run ./0_scripts/c_ring_aln_fun.py

# Reading PDB IDs and their corresponding c-ring chains' IDs from pdb_data file
pdb_data = list(csv.reader(open('pdb_data.csv')))[1:]
pdb_data = sorted(pdb_data, key=lambda pdb: pdb[0])

# Load and clean input structures
cmd.loadall('./1_input_pdb/*')
remove not polymer.protein
###remove resn UNK

# Construct three orthogonal pseudoatom vectors defining the reference XYZ coordinates for ATPase alignment
# 50 points along the positive x and y axes
# 100 points along the postive z axis
for x in range(0, 50, 1): cmd.pseudoatom('X_AXIS', name='PSD', resi=str(x+1), resn='ORG', chain='ZZ', color='red', pos=[x, 0, 0])
for y in range(0, 50, 1): cmd.pseudoatom('Y_AXIS', name='S'+str(y+1), resi=str(y+1), resn='STB', chain='ZZ', color='green', pos=[0, y, 0])
for z in range(0, 100, 1): cmd.pseudoatom('Z_AXIS', name='CRC', resi='0', resn='CRC', chain='ZZ', color='blue', pos=[0, 0, z])

# Draw a line of pseudoatoms at the center of each c-ring
for pdb in pdb_data: c_ring_pseudoatoms(pdb)

# Align the central lines of the c-rings to the (+ve) z-axis
for pdb in pdb_data: c_ring_align(pdb)

# Draw a pseudoatom at the center of each stator's atomic coordinates
# then adjust the pseudoatom's z-coordinate value to match the first atom of the c-ring's central line
for pdb in pdb_data: stator_pseudoatoms(pdb)

# Refine the alignment of the c-rings' central lines to the (+ve) z-axis
# while re-orienting all stators' baselines to match the (+ve) y-axis
for pdb in pdb_data: c_ring_align(pdb)

# Refine the alignment of the F1 domains to 8H9S
## Caution: This will mess up the stator orientations of the mobile structures!!
###extra_fit (name CA or name CRC) , 8H9S, cealign      # Note to Islam: Add this user option to the alignment program later.

# Adjust PyMOL's object panel
cmd.order('*', sort='yes')
cmd.group('Axes', members=f"*_AXIS", action='auto', quiet=1)

# Adjust PyMOL's scene
show nb_spheres, X_AXIS or Y_AXIS or Z_AXIS
show spheres, (resn CRC or resn STB) and (not Axes)
reset
turn x, 90

# Save PyMOL scene (with L-key alignment vectors)
save atpase_alignment-pre_clean_up.pse

# Clean up c-rings
for pdb in pdb_data: c_ring_cleanup(pdb)

# Keep only c-rings 
for pdb in pdb_data: keep_only_c_ring(pdb)

# Save ATPases in separate PDB files
for pdb in pdb_data: cmd.save('./2_aligned_pdb/'+pdb[0]+'_aligned.pdb', pdb[0].lower())

# Save PyMOL scene (clean)
save atpase_alignment.pse
