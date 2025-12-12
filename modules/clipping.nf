/*
 * CLIPPING - Remove soft/hard clipped reads (plasma mode only)
 */

process CLIPPING {
    tag "${meta.sample_id}"
    label 'process_low'
    
    input:
    tuple val(meta), path(bam), path(bai)
    
    output:
    tuple val(meta), path("${meta.sample_id}_sorted_3.bam"), path("${meta.sample_id}_sorted_3.bam.bai"), emit: clipped
    
    script:
    """
    # Filter out reads with soft/hard clipping (CIGAR contains S or H)
    samtools view -h ${bam} | \\
        awk 'BEGIN{OFS="\\t"} /^@/ {print; next} \$6 !~ /[SH]/ {print}' | \\
        samtools view -b -o ${meta.sample_id}_noclip.bam -
    
    samtools sort -@ ${task.cpus} -o ${meta.sample_id}_sorted_3.bam ${meta.sample_id}_noclip.bam
    samtools index ${meta.sample_id}_sorted_3.bam
    """
}
