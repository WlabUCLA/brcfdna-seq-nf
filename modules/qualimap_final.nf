/*
 * QC_FINAL - Final quality control with Qualimap
 */

process QC_FINAL {
    tag "${meta.sample_id}"
    label 'process_medium'
    
    publishDir "${params.outdir}/qc/final/${meta.sample_id}", mode: 'copy'
    
    input:
    tuple val(meta), path(nuc_bam), path(mit_bam)
    
    output:
    tuple val(meta), path("${meta.sample_id}_qc_final/"), emit: qc_dir
    
    script:
    def mem = params.qualimap_mem ?: '15G'
    """
    # Index if needed
    if [[ ! -f ${nuc_bam}.bai ]]; then
        samtools index ${nuc_bam}
    fi
    
    qualimap bamqc \\
        -bam ${nuc_bam} \\
        -outdir ${meta.sample_id}_qc_final \\
        --java-mem-size=${mem} \\
        -nt ${task.cpus}
    """
}
