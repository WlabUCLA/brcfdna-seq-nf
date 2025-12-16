/*
================================================================================
 * BBMERGE
================================================================================
Merge paired-end reads into single cfDNA fragments
 */

process BBMERGE {
    tag "${meta.sample_id}"
    label 'process_medium'
    
    publishDir "${params.outdir}/merged/${meta.sample_id}", mode: 'copy', pattern: '*.log'
    
    input:
    tuple val(meta), path(r1), path(r2)
    
    output:
    tuple val(meta), path("${meta.sample_id}.merged.fastq.gz"), emit: merged
    tuple val(meta), path("${meta.sample_id}.bbmerge.log"), emit: log
    
    script:
    """
    bbmerge.sh \\
        in1=${r1} \\
        in2=${r2} \\
        out=${meta.sample_id}.merged.fastq.gz \\
        adapters=default \\
        mininsert=30 \\
        t=${task.cpus} \\
        2>&1 | tee ${meta.sample_id}.bbmerge.log
    """
}
