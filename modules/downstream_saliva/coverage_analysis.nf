/*
================================================================================
    GENOME COVERAGE ANALYSIS MODULE
================================================================================
    Calculates genome-wide coverage in fixed-size bins using deeptools.
    Useful for copy number analysis and coverage uniformity assessment.
*/

process COVERAGE_ANALYSIS {
    tag "${meta.sample_id}"
    label 'process_medium'
    
    publishDir "${params.outdir}/downstream/coverage/${meta.sample_id}", mode: 'copy'
    
    input:
    tuple val(meta), path(bam)
    
    output:
    tuple val(meta), path("${meta.sample_id}_multibam.csv"), emit: coverage_csv
    tuple val(meta), path("${meta.sample_id}_multibam.npz"), emit: coverage_npz, optional: true
    
    script:
    def bin_size = params.coverage_bin_size ?: 1000000
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    # Index BAM if needed
    if [[ ! -f ${bam}.bai ]]; then
        samtools index ${bam}
    fi
    
    # Run multiBamSummary for coverage binning
    multiBamSummary bins \\
        --bamfiles ${bam} \\
        --binSize ${bin_size} \\
        --minMappingQuality 10 \\
        --outRawCounts ${meta.sample_id}_multibam.csv \\
        -o ${meta.sample_id}_multibam.npz
    
    echo "Coverage analysis complete for ${meta.sample_id}"
    """
}
