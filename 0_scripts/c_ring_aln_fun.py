import csv
from pymol import cmd
from pymol import stored
import numpy as np

def construct_bioassembly(pdb_list, bio_assembly):
    for pdb in pdb_list:
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

def c_ring_pseudoatoms(pdb):
    pdb_id = pdb[0].lower()
    print('Creating a c-ring alignment vector for '+pdb_id)
    cmd.select(pdb_id+'_c_ring', '('+'chain '+pdb[1].replace(',','+')+') and '+pdb_id)
    chain_ids = pdb[1].split(',')
    myspace = {'residues': []}          # This is a user defined PyMOL namespace dictionary
    cmd.iterate_state(-1, pdb_id+' and polymer.protein and chain '+chain_ids[0], 'residues.append(resv)', space=myspace)
    seq_start = sorted(set(myspace['residues']))[0]
    seq_end = sorted(set(myspace['residues']))[-1]
    for res in range(seq_start, (seq_end + 1), 1): cmd.pseudoatom(pdb_id, name='CRC', resi='1', resn='CRC', chain='ZZ', color='hotpink', pos=centroid(pdb_id+'_c_ring and resi '+str(res)+' and name CA'))

def stator_pseudoatoms(pdb):
    pdb_id = pdb[0].lower()
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
        coord1 = cmd.get_coords(f'{pdb_id}_CRC_1st_atom')[0]
        
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

def c_ring_align(pdb):
    cmd.select(pdb[0].lower()+'_c_ring_L_keys', '(resn CRC or resn STB) and '+pdb[0])
    cmd.align(pdb[0].lower()+'_c_ring_L_keys', 'Y_AXIS or Z_AXIS')

def c_ring_cleanup(pdb):
    cmd.remove(pdb[0].lower()+'_c_ring_L_keys')
    cmd.delete(pdb[0].lower()+'_c_ring_L_keys')

def keep_only_c_ring(pdb):
    cmd.remove(pdb[0].lower()+' and not '+pdb[0].lower()+'_c_ring')
    cmd.delete(pdb[0].lower()+'_stator or '+pdb[0].lower()+'_CRC_1st_atom or '+pdb[0].lower()+'_cyto_domain')


cmd.extend("construct_bioassembly", construct_bioassembly)
cmd.extend("c_ring_pseudoatoms", c_ring_pseudoatoms)
cmd.extend("stator_pseudoatoms", stator_pseudoatoms)
cmd.extend("c_ring_align", c_ring_align)
cmd.extend("c_ring_cleanup", c_ring_cleanup)
cmd.extend("keep_only_c_ring", keep_only_c_ring)
