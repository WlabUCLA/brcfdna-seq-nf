/*
================================================================================
 * SAM2BAM
================================================================================
Convert SAM to sorted, indexed BAM
 */

process SAM2BAM {
    tag "${meta.sample_id}"
    label 'process_medium'
    
    input:
    tuple val(meta), path(sam)
    
    output:
    tuple val(meta), path("${meta.sample_id}_sorted_1.bam"), path("${meta.sample_id}_sorted_1.bam.bai"), emit: bam
    
    script:
    """
    samtools sort \\
        -@ ${task.cpus} \\
        -o ${meta.sample_id}_sorted_1.bam \\
        ${sam}
    
    samtools index ${meta.sample_id}_sorted_1.bam
    
    rm -f ${sam}
    """
}
