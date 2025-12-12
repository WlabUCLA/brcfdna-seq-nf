/*
 * ALIGN_BT2 - Bowtie2 alignment for saliva cfDNA
 */

process ALIGN_BT2 {
    tag "${meta.sample_id}"
    label 'process_high'
    
    input:
    tuple val(meta), path(trimmed_fastq)
    
    output:
    tuple val(meta), path("${meta.sample_id}.bowtie2.sam"), emit: sam
    
    script:
    if (!params.bowtie2_index) {
        error "ERROR: --bowtie2_index is required for saliva mode alignment"
    }
    """
    bowtie2 \\
        -x ${params.bowtie2_index} \\
        -U ${trimmed_fastq} \\
        -p ${task.cpus} \\
        --very-sensitive \\
        -S ${meta.sample_id}.bowtie2.sam
    """
}
