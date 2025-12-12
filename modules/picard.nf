/*
 * PICARD - Mark and remove duplicates
 */

process PICARD {
    tag "${meta.sample_id}"
    label 'process_medium'
    
    publishDir "${params.outdir}/deduped/${meta.sample_id}", mode: 'copy', pattern: '*.txt'
    
    input:
    tuple val(meta), path(bx_bam), path(umi_tsv)
    
    output:
    tuple val(meta), path("${meta.sample_id}_sorted_2.bam"), path("${meta.sample_id}_sorted_2.bam.bai"), path("${meta.sample_id}_deduped.metrics.txt"), emit: deduped
    
    script:
    """
    # Sort by coordinate before deduplication (required for MarkDuplicates)
    samtools sort -@ ${task.cpus} -o ${meta.sample_id}_presort.bam ${bx_bam}
    
    picard MarkDuplicates \\
        I=${meta.sample_id}_presort.bam \\
        O=${meta.sample_id}_deduped.bam \\
        M=${meta.sample_id}_deduped.metrics.txt \\
        REMOVE_DUPLICATES=true \\
        ASSUME_SORTED=true
    
    samtools sort -@ ${task.cpus} -o ${meta.sample_id}_sorted_2.bam ${meta.sample_id}_deduped.bam
    samtools index ${meta.sample_id}_sorted_2.bam
    """
}
