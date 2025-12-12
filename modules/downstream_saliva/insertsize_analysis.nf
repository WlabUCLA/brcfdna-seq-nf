/*
================================================================================
    INSERT SIZE DISTRIBUTION ANALYSIS MODULE
================================================================================
    Calculates the distribution of fragment lengths (insert sizes) in cfDNA.
    Fragment size profiles are a key biomarker in liquid biopsy analysis.
*/

process INSERTSIZE_ANALYSIS {
    tag "${meta.sample_id}"
    label 'process_low'
    
    publishDir "${params.outdir}/downstream/insertsize/${meta.sample_id}", mode: 'copy'
    
    input:
    tuple val(meta), path(bam)
    
    output:
    tuple val(meta), path("${meta.sample_id}_insertsize_raw.csv"), emit: insertsize_raw
    tuple val(meta), path("${meta.sample_id}_insertsize_sorted.csv"), emit: insertsize_sorted
    tuple val(meta), path("${meta.sample_id}_insertsize.csv"), emit: insertsize_filled
    
    script:
    def max_size = params.insertsize_max ?: 300
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    # Extract insert sizes from BAM
    # Using samtools to get fragment lengths
    samtools view ${bam} | \\
        awk '{if(\$9 > 0 && \$9 <= ${max_size}) print \$9}' | \\
        sort -n | \\
        uniq -c | \\
        awk 'BEGIN{print "count,size"} {print \$1","\$2}' > ${meta.sample_id}_insertsize_raw.csv
    
    # Parse and process insert size data
    python3 ${projectDir}/bin/insert_size.py \\
        ${meta.sample_id}_insertsize_raw.csv \\
        ${meta.sample_id}_insertsize_sorted.csv \\
        ${meta.sample_id}_insertsize.csv
    
    echo "Insert size analysis complete for ${meta.sample_id}"
    """
}
