/*
================================================================================
 * PICARD - Mark and remove duplicates (UMI-aware)
================================================================================
Uses Picard MarkDuplicates with BX barcode tag for UMI-aware deduplication.
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
    samtools sort -@ ${task.cpus} -o ${meta.sample_id}_sorted_BX.bam ${bx_bam}
    samtools index ${meta.sample_id}_sorted_BX.bam
    
    # UMI-aware duplicate removal using BX tag
    picard MarkDuplicates \
        BARCODE_TAG=BX \
        VALIDATION_STRINGENCY=LENIENT \
        I=${meta.sample_id}_sorted_BX.bam \
        O=${meta.sample_id}_deduped.bam \
        M=${meta.sample_id}_deduped.metrics.txt \
        REMOVE_DUPLICATES=TRUE
    
    # Final sort and index
    samtools sort -@ ${task.cpus} -o ${meta.sample_id}_sorted_2.bam ${meta.sample_id}_deduped.bam
    samtools index ${meta.sample_id}_sorted_2.bam
    
    # Clean up intermediate files
    rm -f ${meta.sample_id}_sorted_BX.bam ${meta.sample_id}_sorted_BX.bam.bai ${meta.sample_id}_deduped.bam
    """
}
