/*
================================================================================
    END MOTIF ANALYSIS MODULE
================================================================================
    Extracts and analyzes the first 4 base pairs (end motifs) from each read.
    Useful for nuclease signature detection and fragmentation pattern analysis.
*/

process ENDMOTIF_ANALYSIS {
    tag "${meta.sample_id}"
    label 'process_medium'
    
    publishDir "${params.outdir}/downstream/endmotif/${meta.sample_id}", mode: 'copy'
    
    input:
    tuple val(meta), path(bam)
    
    output:
    tuple val(meta), path("${meta.sample_id}_endmotif.csv"), emit: motif_csv
    tuple val(meta), path("${meta.sample_id}_endmotif_ordered.csv"), emit: motif_ordered
    tuple val(meta), path("${meta.sample_id}_first4.csv"), emit: first4_csv, optional: true
    
    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    # Index BAM if needed
    if [[ ! -f ${bam}.bai ]]; then
        samtools index ${bam}
    fi
    
    # Extract first 4 base pairs and count motifs
    python3 ${projectDir}/bin/motif_count.py \\
        ${bam} \\
        ${meta.sample_id}_first4.csv \\
        ${meta.sample_id}_endmotif.csv
    
    # Order motifs alphabetically and fill missing combinations
    python3 ${projectDir}/bin/motif_ordered.py \\
        ${meta.sample_id}_endmotif.csv \\
        ${meta.sample_id}_endmotif_ordered.csv
    
    echo "End motif analysis complete for ${meta.sample_id}"
    """
}
