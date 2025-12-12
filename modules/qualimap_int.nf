/*
 * QC_INT - Intermediate quality control with Qualimap
 */

process QC_INT {
    tag "${meta.sample_id}"
    label 'process_medium'
    
    publishDir "${params.outdir}/qc/intermediate/${meta.sample_id}", mode: 'copy'
    
    input:
    tuple val(meta), path(bam), path(bai)
    
    output:
    tuple val(meta), path("${meta.sample_id}_qc_int/"), emit: qc_dir
    
    script:
    def mem = params.qualimap_mem ?: '15G'
    """
    qualimap bamqc \\
        -bam ${bam} \\
        -outdir ${meta.sample_id}_qc_int \\
        --java-mem-size=${mem} \\
        -nt ${task.cpus}
    """
}
