#!/usr/bin/env nextflow

process pymol_c_ring_aligner {

    publishDir "${projectDir}/results/atpase/${cif.baseName}", mode: 'copy', pattern: '*_aligned_atpase.pdb'
    publishDir "${projectDir}/results/rotor/${cif.baseName}", mode: 'copy', pattern: '*_aligned_rotor.pdb'
    publishDir "${projectDir}/results/atpase/${cif.baseName}", mode: 'copy', pattern: '*_rotor_data.csv'
    publishDir "${projectDir}/results/rotor/${cif.baseName}", mode: 'copy', pattern: '*_rotor_data.csv'

    input:
    val pymol_c_ring_aligner
    val scripts_dir
    path cif

    output:
    path '*_aligned_atpase.pdb' , emit: aligned_atpase
    path '*_aligned_rotor.pdb' , emit: aligned_rotor
    path '*_rotor_data.csv'

    script:
    """
    python ${pymol_c_ring_aligner} ${scripts_dir} ${cif}
    """
}

process atpase_pdbfixer {

    label 'pdbfixer'
    publishDir "${projectDir}/results/atpase/${aligned_atpase.baseName.substring(0, 4)}", mode: 'copy'

    input:
    path aligned_atpase
    val platform

    output:
    path '*_fixed_atpase.pdb' , emit: fixed_atpase

    script:
    """
    export OPENMM_DEFAULT_PLATFORM=${platform}
    echo "Fixing ${aligned_atpase.baseName}"
    pdbfixer "${aligned_atpase}" \
        --output="${aligned_atpase.baseName.substring(0, 4)}_fixed_atpase.pdb" \
        --add-atoms=heavy \
        --keep-heterogens=all \
        --replace-nonstandard \
        --verbose
    """
}

process rotor_pdbfixer {

    label 'pdbfixer'
    publishDir "${projectDir}/results/rotor/${aligned_rotor.baseName.substring(0, 4)}", mode: 'copy'

    input:
    path aligned_rotor
    val platform

    output:
    path '*_fixed_rotor.pdb' , emit: fixed_rotor

    script:
    """
    export OPENMM_DEFAULT_PLATFORM=${platform}
    echo "Fixing ${aligned_rotor.baseName}"
    pdbfixer "${aligned_rotor}" \
        --output="${aligned_rotor.baseName.substring(0, 4)}_fixed_rotor.pdb" \
        --add-atoms=heavy \
        --keep-heterogens=all \
        --replace-nonstandard \
        --verbose
    """
}

process atpase_pdb2pqr {

    label 'pdb2pqr'
    publishDir "${projectDir}/results/atpase/${fixed_atpase.baseName.substring(0, 4)}", mode: 'copy'

    input:
    val pdb2pqr_script
    val scripts_dir
    path fixed_atpase
    val pH

    output:
    path "*.pqr" , emit: atpase_pqr
    path "*.in" , emit: atpase_apbs_in
    path "*_pdb2pqr.log"

    script:
    """
    python ${pdb2pqr_script} ${scripts_dir} ${fixed_atpase} ${pH}
    """
}

process rotor_pdb2pqr {

    label 'pdb2pqr'
    publishDir "${projectDir}/results/rotor/${fixed_rotor.baseName.substring(0, 4)}", mode: 'copy'

    input:
    val pdb2pqr_script
    val scripts_dir
    path fixed_rotor
    val pH

    output:
    path "*.pqr" , emit: rotor_pqr
    path "*.in" , emit: rotor_apbs_in
    path "*_pdb2pqr.log"

    script:
    """
    python ${pdb2pqr_script} ${scripts_dir} ${fixed_rotor} ${pH}
    """
}

process atpase_ct_scan {

    label 'apbs'
    publishDir "${projectDir}/results/atpase/${atpase_pqr.baseName.substring(0, 4)}", mode: 'copy'

    input:
    path atpase_pqr
    path atpase_apbs_in
    val apbs_ctscan_analyzer
    val distmin
    val distmax

    output:
    path "*.dx" , emit: atpase_mesp
    path "*.mc" , emit: atpase_mc
    path "*_apbs.pqr"
    path "*.dat"
    path "pot_avg"
    path "pot_ct_scan"

    script:
    """
    apbs ${atpase_apbs_in}
    mv io.mc "${atpase_pqr.baseName}_io.mc"
    python ${apbs_ctscan_analyzer} \\
        --pqr ${atpase_pqr} \\
        --dx ${atpase_pqr.baseName}-pot.dx \\
        --distmin ${distmin} \\
        --distmax ${distmax}
    """
}

process rotor_ct_scan {

    label 'apbs'
    publishDir "${projectDir}/results/rotor/${rotor_pqr.baseName.substring(0, 4)}", mode: 'copy'

    input:
    path rotor_pqr
    path rotor_apbs_in
    val apbs_ctscan_analyzer
    val distmin
    val distmax

    output:
    path "*.dx" , emit: rotor_mesp
    path "*.mc" , emit: rotor_mc
    path "*_apbs.pqr"
    path "*.dat"
    path "pot_avg"
    path "pot_ct_scan"

    script:
    """
    apbs ${rotor_apbs_in}
    mv io.mc "${rotor_pqr.baseName}_io.mc"
    python ${apbs_ctscan_analyzer} \\
        --pqr ${rotor_pqr} \\
        --dx ${rotor_pqr.baseName}-pot.dx \\
        --distmin ${distmin} \\
        --distmax ${distmax}
    """
}
