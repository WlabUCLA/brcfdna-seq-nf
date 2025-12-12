/*
================================================================================
    NUCLEASE SIGNATURE ANALYSIS MODULE
================================================================================
    Runs comprehensive nuclease signature analysis including:
    - Nuclease cleavage patterns
    - DNA shape analysis
    - Structure detection
    - Damage analysis
*/

process NUCLEASE_SIGNATURE_ANALYSIS {
    tag "${meta.sample_id}"
    label 'process_high'
    
    publishDir "${params.outdir}/downstream/nuclease_signature/${meta.sample_id}", mode: 'copy'
    
    input:
    tuple val(meta), path(bam)
    path reference
    
    output:
    tuple val(meta), path("${meta.sample_id}_nuclease_signature/"), emit: nuclease_results
    tuple val(meta), path("${meta.sample_id}_nuclease_signature/comprehensive_analysis_results.json"), emit: nuclease_json, optional: true
    tuple val(meta), path("${meta.sample_id}_nuclease_signature/nuclease_analysis_comprehensive.csv"), emit: nuclease_csv, optional: true
    tuple val(meta), path("${meta.sample_id}_nuclease_signature/comprehensive_analysis_summary.csv"), emit: summary_csv, optional: true
    
    script:
    def ref_arg = reference.size() > 0 ? "-r ${reference}" : ""
    def fast_arg = params.nuclease_fast_mode ? "--fast-mode" : ""
    def chunk_arg = params.nuclease_chunk_size ? "--chunk-size ${params.nuclease_chunk_size}" : ""
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    # Index BAM if needed
    if [[ ! -f ${bam}.bai ]]; then
        samtools index ${bam}
    fi
    
    # Run nuclease signature analysis
    python3 ${projectDir}/bin/nuclease_signature.py \\
        -i ${bam} \\
        ${ref_arg} \\
        -o ${meta.sample_id}_nuclease_signature \\
        ${fast_arg} \\
        ${chunk_arg}
    
    echo "Nuclease signature analysis complete for ${meta.sample_id}"
    """
}
