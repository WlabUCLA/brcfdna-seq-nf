/*
================================================================================
 * FASTP_TRIM
================================================================================
Adapter trimming and quality filtering
 */

process FASTP_TRIM {
    tag "${meta.sample_id}"
    label 'process_medium'
    
    publishDir "${params.outdir}/trimmed/${meta.sample_id}", mode: 'copy', pattern: '*.{html,json}'
    
    input:
    tuple val(meta), path(merged_fastq)
    
    output:
    tuple val(meta), path("${meta.sample_id}.trimmed.fastq.gz"), path("${meta.sample_id}.fastp.html"), path("${meta.sample_id}.fastp.json"), emit: trimmed
    
    script:
    """
    fastp \\
        -i ${merged_fastq} \\
        -o ${meta.sample_id}.trimmed.fastq.gz \\
        --html ${meta.sample_id}.fastp.html \\
        --json ${meta.sample_id}.fastp.json \\
        --thread ${task.cpus} \\
        --length_required 30 \\
        --qualified_quality_phred 20
    """
}
