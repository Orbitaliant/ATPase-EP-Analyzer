#!/usr/bin/env nextflow

// Pipeline modules
include {pymol_c_ring_aligner} from "${projectDir}/0_scripts/atpase_ep_analyzer_proc.nf"
include {atpase_pdbfixer; rotor_pdbfixer} from "${projectDir}/0_scripts/atpase_ep_analyzer_proc.nf"
include {atpase_pdb2pqr; rotor_pdb2pqr} from "${projectDir}/0_scripts/atpase_ep_analyzer_proc.nf"
include {atpase_ct_scan; rotor_ct_scan} from "${projectDir}/0_scripts/atpase_ep_analyzer_proc.nf"

// CIF input directory
params.cif_input_dir = "${projectDir}/cif_input"

// Workflow
workflow {
    cif_ch = channel.fromPath("${params.cif_input_dir}/*.cif")

    pymol_c_ring_aligner(params.pymol_c_ring_aligner, params.scripts_dir, cif_ch)

    atpase_pdbfixer(pymol_c_ring_aligner.out.aligned_atpase, params.openmm_platform)

    rotor_pdbfixer(pymol_c_ring_aligner.out.aligned_rotor, params.openmm_platform)

    atpase_pdb2pqr(params.pdb2pqr_script, params.scripts_dir, atpase_pdbfixer.out.fixed_atpase, params.pH)

    rotor_pdb2pqr(params.pdb2pqr_script, params.scripts_dir, rotor_pdbfixer.out.fixed_rotor, params.pH)

    atpase_ct_scan(atpase_pdb2pqr.out.atpase_pqr, atpase_pdb2pqr.out.atpase_apbs_in, params.apbs_ctscan_analyzer, params.distmin, params.distmax)

    rotor_ct_scan(rotor_pdb2pqr.out.rotor_pqr, rotor_pdb2pqr.out.rotor_apbs_in, params.apbs_ctscan_analyzer, params.distmin, params.distmax)
}
