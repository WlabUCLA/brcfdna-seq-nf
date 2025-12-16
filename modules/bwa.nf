/*
================================================================================
 * ALIGN_BWA
================================================================================
BWA-MEM alignment for plasma cfDNA
 */

process ALIGN_BWA {
    tag "${meta.sample_id}"
    label 'process_high'
    
    publishDir "${params.outdir}/aligned/${meta.sample_id}", mode: 'copy', pattern: '*.bai'
    
    input:
    tuple val(meta), path(trimmed_fastq)
    
    output:
    tuple val(meta), path("${meta.sample_id}_sorted_1.bam"), path("${meta.sample_id}_sorted_1.bam.bai"), emit: bam
    
    script:
    if (!params.bwa_index) {
        error "ERROR: --bwa_index is required for plasma mode alignment"
    }
    """
    bwa mem \\
        -t ${task.cpus} \\
        -M \\
        ${params.bwa_index} \\
        ${trimmed_fastq} | \\
    samtools sort -@ ${task.cpus} -o ${meta.sample_id}_sorted_1.bam -
    
    samtools index ${meta.sample_id}_sorted_1.bam
    """
}
