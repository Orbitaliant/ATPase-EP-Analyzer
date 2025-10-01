import sys
import csv
from mmcif.io.IoAdapterCore import IoAdapterCore
from pymol import cmd
from pymol import stored
from centroid import centroid
import numpy as np
import pymol

# ==============================
# FUNCTION: MMCIF Rotor Ring Detector
# ==============================

io = IoAdapterCore()

def rotor_ring_detector(cif):
    # Extract PDB ID from cif path
    pdb_id = cif.split("/")[-1][:-4].lower()

    # Extract rotor helices chain identifiers
    list_data_container = io.readFile(cif)
    data_container = list_data_container[0]
    c_ring_helices_ids = None  # Default value
    for i in range(len(data_container.getObj('entity_poly'))):
        entity = data_container.getObj('entity')[i]
        entity_poly = data_container.getObj('entity_poly')[i]
        if int(entity[5]) >= 8:
            c_ring_helices_ids = entity_poly[6]

    # Write sorted output to CSV
    output_csv = f'{pdb_id}_rotor_data.csv'
    with open(output_csv, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow([pdb_id, c_ring_helices_ids])

    print(f"{output_csv} created successfully.")

# ==============================
# FUNCTION: PyMOL Construct Cartesian Coordinates Using Pseudoatoms
# ==============================
def construct_cartesian_coordinates_pseudoatoms(xlen=50, ylen=50, zlen=100):
    for x in range(0, xlen, 1): cmd.pseudoatom('X_AXIS', name='PSD', resi=str(x+1), resn='ORG', chain='ZZ', color='red', pos=[x, 0, 0])
    for y in range(0, ylen, 1): cmd.pseudoatom('Y_AXIS', name='S'+str(y+1), resi=str(y+1), resn='STB', chain='ZZ', color='green', pos=[0, y, 0])
    for z in range(0, zlen, 1): cmd.pseudoatom('Z_AXIS', name='CRC', resi='0', resn='CRC', chain='ZZ', color='blue', pos=[0, 0, z])

cmd.extend("construct_cartesian_coordinates_pseudoatoms", construct_cartesian_coordinates_pseudoatoms)

# ==============================
# FUNCTION: PyMOL Construct C-Ring Bioassembly (optional)
# ==============================
# only use if the assymetric unit of the pdb entry containing the c-ring consists of only one c-subunit
# relies on symmetry to reconstruct the c-ring
def construct_bioassembly(pdb, bio_assembly):
    cmd.set('assembly', f'{str(bio_assembly)}')
    cmd.fetch(pdb, type='cif')
    cmd.split_states(pdb)
    cmd.delete(pdb)
    assembly_parts = cmd.get_names('all')
    for i in range(1, len(assembly_parts)+1):
        cmd.alter(assembly_parts[i-1], f'chain=chain+str({i})')
    cmd.create(pdb, 'all')
    cmd.delete(f'not {pdb}')
    cmd.save(f'{pdb}.cif', format='cif')

cmd.extend("construct_bioassembly", construct_bioassembly)

# ==============================
# FUNCTION: PyMOL Construct C-Ring Central Axis Using Pseudoatoms
# ==============================
def c_ring_pseudoatoms(pdb):
    pdb_id = pdb[0]
    print(f'Creating a c-ring alignment vector for {pdb_id}')
    cmd.select(f"{pdb_id}_c_ring", f"(chain {pdb[1].replace(',', '+')}) and {pdb_id}")
    chain_ids = pdb[1].split(',')
    myspace = {'residues': []}
    cmd.iterate_state(-1, f"{pdb_id} and polymer.protein and chain {chain_ids[0]}", 'residues.append(resv)', space=myspace)
    seq_start = sorted(set(myspace['residues']))[0]
    seq_end = sorted(set(myspace['residues']))[-1]
    for res in range(seq_start, seq_end + 1):
        sel = f"{pdb_id}_c_ring and resi {res} and name CA"
        n = cmd.count_atoms(sel)
        if n == 0:
            print(f"WARNING: No CA atoms found for resi {res} ({sel}). Skipping.")
            continue
        try:
            pos = centroid(sel)
        except Exception as e:
            print(f"ERROR: Could not compute centroid for {sel}: {e}")
            continue
        cmd.pseudoatom(pdb_id, name='CRC', resi='1', resn='CRC', chain='ZZ', color='hotpink', pos=pos)

cmd.extend("c_ring_pseudoatoms", c_ring_pseudoatoms)

# ==============================
# FUNCTION: PyMOL Construct C-Ring-Stator-Base Radial Vector Using Pseudoatoms
# ==============================
def stator_pseudoatoms(pdb):
    pdb_id = pdb[0]
    print('Creating a stator alignment vector for ' + pdb_id)

    # Select cytoplasmic domain
    cmd.select(pdb_id + '_cyto_domain', f'(bychain {pdb_id} near_to 75 of {pdb_id}_c_ring) and ((x > -30 and x < 30) and (y > -30 and y < 30) and (z > 50))')

    # Refine selection to include all atoms in chains near the cytoplasmic domain
    cmd.select(pdb_id + '_cyto_domain', f'bychain {pdb_id}_cyto_domain')

    # Select the stator, excluding the c-ring and cytoplasmic domain
    cmd.select(pdb_id + '_stator', f'({pdb_id} and (not {pdb_id}_c_ring) and (not {pdb_id}_cyto_domain))')

    # Check for ATPases with missing stators
    if cmd.count_atoms(pdb_id + '_stator') == 0:
        print(f'\nWarning! PDB structure {pdb_id} has no stator domain...\n')
    
    else:
    
        # Create a pseudoatom at the centroid of the stator selection
        cmd.pseudoatom(pdb_id, name='PSD', resn='TMP', chain='XY', color='hotpink', pos=centroid(f'{pdb_id}_stator and name CA'))

        # Get the coordinates of the first atom of the c-ring central line
        cmd.select(pdb_id + '_CRC_1st_atom', 'first (resn CRC) and '+pdb_id)
        coord1 = cmd.get_coords(f'{pdb_id}_CRC_1st_atom')
        if coord1 is None:
            print(f"Warning: No coordinates found for {pdb_id}_CRC_1st_atom")
            return
        coord1 = coord1[0]
        
        # Align the pseudoatom along the z-axis
        cmd.alter_state(1, f'{pdb_id} and chain XY and resn TMP', f'z = {coord1[2]}')
        
        # Rebuild the view to reflect changes
        cmd.rebuild(representation='everything')
    
        # Get the coordinates of the stator's z-aligned centroid
        coord2 = cmd.get_coords(f'{pdb_id} and resn TMP')[0]
    
        # Calculate the vector between the coordinates of the two atoms
        vector = np.array(coord2) - np.array(coord1)
    
        # Normalize the vector to 1 Ångström length
        vector_normalized = vector / np.linalg.norm(vector)
    
        # Number of dummy atoms to create
        num_dummy_atoms = 50
    
        # Generate coordinates for the dummy atoms along the vector
        dummy_coords = [coord1 + i * vector_normalized for i in range(num_dummy_atoms)]
    
        # Create dummy atoms at each generated coordinate
        exec("; ".join([f"cmd.pseudoatom(object='{pdb_id}', name='S{i}', resi='{i}', resn='STB', chain='ZZ', color='yellow', pos={list(coord)})" for i, coord in enumerate(dummy_coords, 1)]))

        # Remove all stators' temporary centroid pseudoatoms
        cmd.remove("resn TMP")

cmd.extend("stator_pseudoatoms", stator_pseudoatoms)

# ==============================
# FUNCTION: PyMOL Align C-Ring Pseudoatom L-Keys To Pseudoatom Cartesian Coordinates
# ==============================
def c_ring_align(pdb):
    cmd.select(pdb[0]+'_c_ring_L_keys', '(resn CRC or resn STB) and '+pdb[0])
    cmd.align(pdb[0]+'_c_ring_L_keys', 'Y_AXIS or Z_AXIS')

cmd.extend("c_ring_align", c_ring_align)

# ==============================
# FUNCTION: PyMOL Delete C-Ring Pseudoatom L-Keys
# ==============================
def c_ring_cleanup(pdb):
    cmd.remove(pdb[0]+'_c_ring_L_keys')
    cmd.delete(pdb[0]+'_c_ring_L_keys')

cmd.extend("c_ring_cleanup", c_ring_cleanup)

# ==============================
# FUNCTION: PyMOL Extract Aligned C-Ring
# ==============================
def extract_rotor(pdb):
    cmd.remove(pdb[0]+' and not '+pdb[0]+'_c_ring')
    cmd.delete(pdb[0]+'_stator or '+pdb[0]+'_CRC_1st_atom or '+pdb[0]+'_cyto_domain')

cmd.extend("extract_rotor", extract_rotor)

# ==============================
# FUNCTION: PyMOL C-Ring Aligner
# ==============================
def pymol_c_ring_aligner(cif):
    # Start PyMOL in command-line mode
    pymol.finish_launching(['pymol', '-c'])

    # Read rotor data from input csv file
    pdb_id = cif.split("/")[-1][:-4].lower()
    rot_dat = list(csv.reader(open(f'{pdb_id}_rotor_data.csv')))[0]#[1]

    # Load input cif file
    try:
        cmd.load(cif)
    except pymol.CmdException as e:
        sys.exit(f'Failed to load {rot_dat[0]}: {e}')

    # Clean the loaded cif structure
    cmd.remove("not polymer.protein")
    cmd.remove("not alt ''+A")
    cmd.alter("all", "alt=''")

    # Execute c-ring alignment functions
    construct_cartesian_coordinates_pseudoatoms()
    c_ring_pseudoatoms(rot_dat)
    c_ring_align(rot_dat)
    stator_pseudoatoms(rot_dat)
    c_ring_align(rot_dat)
    c_ring_cleanup(rot_dat)

    # Save the aligned atpase
    cmd.save(f'{rot_dat[0]}_aligned_atpase.pdb', rot_dat[0])

    # Extract & save the aligned c-ring
    extract_rotor(rot_dat)
    cmd.save(f'{rot_dat[0]}_aligned_rotor.pdb', rot_dat[0])

    # Ensure PyMOL quits after processing
    cmd.quit()

    print(f'{rot_dat[0]} alignment completed successfully...')

# ============================== #
#               Main             #
# ============================== #

if __name__ == "__main__":
    scripts_dir = sys.argv[1]
    cif = sys.argv[2]
    
    sys.path.append(scripts_dir)
    rotor_ring_detector(cif)
    pymol_c_ring_aligner(cif)

# Dev comments
# - Define pdb_id = rot_dat[0] in pymol_c_ring_aligner
# - Remove pdb_id = pdb[0] from functions
# - Replace pdb[0] with pdb_id in all functions except pymol_c_ring_aligner
# - Define rotor_ids = rot_dat[1] in pymol_c_ring_aligner
# - Replace pdb[1] with rotor_ids in all functions except pymol_c_ring_aligner
