/*
================================================================================
 * BLACKLIST
================================================================================
Remove reads overlapping blacklist regions
 */

process BLACKLIST {
    tag "${meta.sample_id}"
    label 'process_low'
    
    publishDir "${params.outdir}/final/${meta.sample_id}", mode: 'copy'
    
    input:
    tuple val(meta), path(nuc_bam), path(mit_bam)
    
    output:
    tuple val(meta), path("${meta.sample_id}_sorted_blacklisted.bam"), path("${meta.sample_id}_sorted_blacklisted.bam.bai"), emit: blacklisted
    
    script:
    """
    if [[ -n "${params.blacklist_bed}" && -f "${params.blacklist_bed}" ]]; then
        # Remove blacklist regions
        bedtools intersect \\
            -a ${nuc_bam} \\
            -b ${params.blacklist_bed} \\
            -v \\
            > ${meta.sample_id}_blacklisted_unsorted.bam
        
        samtools sort -@ ${task.cpus} -o ${meta.sample_id}_sorted_blacklisted.bam ${meta.sample_id}_blacklisted_unsorted.bam
    else
        # No blacklist - just copy/rename
        cp ${nuc_bam} ${meta.sample_id}_sorted_blacklisted.bam
    fi
    
    samtools index ${meta.sample_id}_sorted_blacklisted.bam
    """
}
