/*
 * UMITOOLS - UMI-based deduplication grouping
 */

process UMITOOLS {
    tag "${meta.sample_id}"
    label 'process_medium'
    
    publishDir "${params.outdir}/umi/${meta.sample_id}", mode: 'copy', pattern: '*.tsv'
    
    input:
    tuple val(meta), path(rx_bam)
    
    output:
    tuple val(meta), path("${meta.sample_id}_BX.bam"), path("${meta.sample_id}_umicorrection.tsv"), emit: bx_bam
    
    script:
    """
    # Sort by position first
    samtools sort -o ${meta.sample_id}_RX_sorted.bam ${rx_bam}
    samtools index ${meta.sample_id}_RX_sorted.bam
    
    umi_tools group \\
        -I ${meta.sample_id}_RX_sorted.bam \\
        --output-bam \\
        -S ${meta.sample_id}_BX.bam \\
        --group-out=${meta.sample_id}_umicorrection.tsv \\
        --extract-umi-method=tag \\
        --umi-tag=RX
    
    rm -f ${meta.sample_id}_RX_sorted.bam*
    """
}
