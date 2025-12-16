/*
================================================================================
    UMITOOLS - UMI-based deduplication grouping
================================================================================
    Groups reads by UMI to identify PCR duplicates.
    Uses directional method for UMI network clustering.
    
    INPUT:  tuple(meta, rx_bam) - BAM with RX tags
    OUTPUT: tuple(meta, bx_bam, tsv) - BAM with BX tags + correction table
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
    # Sort by position first (required for umi_tools)
    samtools sort -o ${meta.sample_id}_sorted_RX.bam ${rx_bam}
    samtools index ${meta.sample_id}_sorted_RX.bam
    
    # UMI grouping with directional method
    umi_tools group \\
        --output-bam \\
        --stdin=${meta.sample_id}_sorted_RX.bam \\
        --stdout=${meta.sample_id}_BX.bam \\
        --group-out=${meta.sample_id}_umicorrection.tsv \\
        --extract-umi-method=tag \\
        --umi-tag=RX \\
        --method=directional \\
        --paired \\
        --unmapped-reads=use \\
        --chimeric-pairs=discard
    
    # Clean up intermediate files
    rm -f ${meta.sample_id}_sorted_RX.bam ${meta.sample_id}_sorted_RX.bam.bai
    """
}
