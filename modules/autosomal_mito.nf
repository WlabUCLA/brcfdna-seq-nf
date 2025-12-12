/*
 * AUTOSOMAL_MITO - Split BAM into nuclear and mitochondrial fractions
 * Accepts input from either CLIPPING (plasma) or PICARD (saliva)
 */

process AUTOSOMAL_MITO {
    tag "${meta.sample_id}"
    label 'process_low'
    
    publishDir "${params.outdir}/split/${meta.sample_id}", mode: 'copy', pattern: '*_mit.bam'
    
    input:
    tuple val(meta), path(bam), path(extra_files)  // extra_files can be bai or [bai, metrics]
    
    output:
    tuple val(meta), path("${meta.sample_id}_nuc.bam"), path("${meta.sample_id}_mit.bam"), emit: split_bams
    
    script:
    """
    # Index if not present
    if [[ ! -f ${bam}.bai ]]; then
        samtools index ${bam}
    fi
    
    # Extract nuclear (autosomal) reads - exclude chrM/MT
    samtools view -h ${bam} | \\
        awk '/^@/ && !/chrM|MT/ {print; next} !/^@/ && \$3 !~ /chrM|MT/ {print}' | \\
        samtools view -b -o ${meta.sample_id}_nuc_unsorted.bam -
    
    # Extract mitochondrial reads
    samtools view -b ${bam} chrM -o ${meta.sample_id}_mit.bam 2>/dev/null || \\
        samtools view -b ${bam} MT -o ${meta.sample_id}_mit.bam 2>/dev/null || \\
        touch ${meta.sample_id}_mit.bam
    
    # Sort nuclear BAM
    samtools sort -@ ${task.cpus} -o ${meta.sample_id}_nuc.bam ${meta.sample_id}_nuc_unsorted.bam
    samtools index ${meta.sample_id}_nuc.bam
    
    rm -f ${meta.sample_id}_nuc_unsorted.bam
    """
}
