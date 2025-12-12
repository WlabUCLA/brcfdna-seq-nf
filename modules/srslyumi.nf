/*
 * SRSLYUMI - Tag reads with UMI (RX tag)
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
    srslyumi tag \\
        -i ${bam} \\
        -o ${meta.sample_id}_RX.bam
    """
}
