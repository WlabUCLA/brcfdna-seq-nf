/*
================================================================================
    ALIGN_BT2 - Bowtie2 alignment for saliva cfDNA
================================================================================
    Aligns trimmed reads to reference genome and separates mapped/unmapped.
    Unmapped reads are saved for downstream microbiome analysis (MetaPhlAn).
    
    OUTPUT:
    - SAM file with mapped reads
    - BAM file with unmapped reads (for MetaPhlAn)
*/

process ALIGN_BT2 {
    tag "${meta.sample_id}"
    label 'process_high'
    
    publishDir "${params.outdir}/unmapped/${meta.sample_id}", mode: 'copy', pattern: '*_unmapped.bam*'
    
    input:
    tuple val(meta), path(trimmed_fastq)
    
    output:
    tuple val(meta), path("${meta.sample_id}.bowtie2.sam"), emit: sam
    tuple val(meta), path("${meta.sample_id}_unmapped.bam"), path("${meta.sample_id}_unmapped.bam.bai"), emit: unmapped
    
    script:
    def bt2_index = params.bowtie2_index
    if (!bt2_index) {
        error "Bowtie2 index not specified. Please provide --bowtie2_index parameter."
    }
    """
    # Align with bowtie2
    bowtie2 \
        -x ${bt2_index} \
        -U ${trimmed_fastq} \
        -p ${task.cpus} \
        --very-sensitive \
        -S ${meta.sample_id}_all.sam
    
    # Extract unmapped reads (flag 4) and save as BAM for MetaPhlAn
    samtools view -f 4 -b ${meta.sample_id}_all.sam | samtools sort -o ${meta.sample_id}_unmapped.bam -
    samtools index ${meta.sample_id}_unmapped.bam
    
    # Extract mapped reads only for downstream processing
    samtools view -F 4 -h ${meta.sample_id}_all.sam > ${meta.sample_id}.bowtie2.sam
    
    # Clean up
    rm -f ${meta.sample_id}_all.sam
    
    # Report stats
    TOTAL=$(samtools view -c ${meta.sample_id}.bowtie2.sam || echo "0")
    UNMAPPED=$(samtools view -c ${meta.sample_id}_unmapped.bam || echo "0")
    echo "Mapped reads: ${TOTAL}"
    echo "Unmapped reads: ${UNMAPPED}"
    """
}
