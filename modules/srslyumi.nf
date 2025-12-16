/*
================================================================================
 * SRSLYUMI
================================================================================
    Uses srslyumi-bamtag to extract UMI from read names and add RX tag.
*/

process SRSLYUMI {
    tag "${meta.sample_id}"
    label 'process_medium'
    
    input:
    tuple val(meta), path(bam), path(bai)
    
    output:
    tuple val(meta), path("${meta.sample_id}_RX.bam"), emit: rx_bam
    
    script:
    """
    srslyumi-bamtag \\
        --binary \\
        -o ${meta.sample_id}_RX.bam \\
        --take-fragment 0 \\
        ${bam}
    """
}
