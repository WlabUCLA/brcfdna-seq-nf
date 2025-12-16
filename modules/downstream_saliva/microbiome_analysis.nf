/*
================================================================================
    MICROBIOME ANALYSIS MODULE (MetaPhlAn 4)
================================================================================
    Performs microbiome profiling on unmapped reads using MetaPhlAn 4.
    Identifies bacterial, archaeal, viral, and eukaryotic species in saliva samples.
    
    INPUT: Unmapped BAM from alignment step (not the final processed BAM)
    
    IMPORTANT: MetaPhlAn requires a pre-installed database.
    
    Database Setup (run once before using this module):
      metaphlan --install --bowtie2db /path/to/metaphlan_db
    
    Then provide the path via:
      --metaphlan_db /path/to/metaphlan_db
*/

process MICROBIOME_ANALYSIS {
    tag "${meta.sample_id}"
    label 'process_high'
    
    publishDir "${params.outdir}/downstream/microbiome/${meta.sample_id}", mode: 'copy'
    
    input:
    tuple val(meta), path(unmapped_bam), path(unmapped_bai)
    
    output:
    tuple val(meta), path("${meta.sample_id}_mpa.csv"), emit: mpa_profile
    tuple val(meta), path("${meta.sample_id}_mpa_bowtie2.bz2"), emit: mpa_bowtie2, optional: true
    tuple val(meta), path("${meta.sample_id}_unmapped.fastq.gz"), emit: unmapped_fastq, optional: true
    
    script:
    def db_arg = params.metaphlan_db ? "--bowtie2db ${params.metaphlan_db}" : ""
    def index_arg = params.metaphlan_index ? "--index ${params.metaphlan_index}" : ""
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    echo "=== MetaPhlAn Microbiome Analysis ==="
    echo "Sample: ${meta.sample_id}"
    echo "Input: ${unmapped_bam}"
    echo "Database: ${params.metaphlan_db ?: 'default (auto-detect)'}"
    
    # Convert unmapped BAM to FASTQ
    echo "Converting unmapped BAM to FASTQ..."
    samtools fastq -@ ${task.cpus} ${unmapped_bam} | gzip > ${meta.sample_id}_unmapped.fastq.gz
    
    # Count unmapped reads
    UNMAPPED_COUNT=\$(zcat ${meta.sample_id}_unmapped.fastq.gz | wc -l)
    UNMAPPED_COUNT=\$((UNMAPPED_COUNT / 4))
    echo "Unmapped reads: \${UNMAPPED_COUNT}"
    
    # Check if we have unmapped reads to process
    if [[ \${UNMAPPED_COUNT} -gt 0 ]]; then
        echo "Running MetaPhlAn 4..."
        
        # Run MetaPhlAn with proper options
        metaphlan ${meta.sample_id}_unmapped.fastq.gz \\
            --input_type fastq \\
            ${db_arg} \\
            ${index_arg} \\
            --nproc ${task.cpus} \\
            --bowtie2out ${meta.sample_id}_mpa_bowtie2.bz2 \\
            --unclassified_estimation \\
            -o ${meta.sample_id}_mpa_raw.txt
        
        # Convert to CSV format
        echo "Converting to CSV format..."
        head -n 1 ${meta.sample_id}_mpa_raw.txt | sed 's/\\t/,/g' > ${meta.sample_id}_mpa.csv
        tail -n +2 ${meta.sample_id}_mpa_raw.txt | sed 's/\\t/,/g' >> ${meta.sample_id}_mpa.csv
        rm -f ${meta.sample_id}_mpa_raw.txt
        
        echo "MetaPhlAn analysis complete"
        
        # Add summary stats
        echo ""
        echo "=== Top Taxa ==="
        head -20 ${meta.sample_id}_mpa.csv
        
    else
        echo "WARNING: No unmapped reads found in ${unmapped_bam}"
        echo "clade_name,relative_abundance" > ${meta.sample_id}_mpa.csv
        echo "UNKNOWN,100.0" >> ${meta.sample_id}_mpa.csv
    fi
    
    echo "Microbiome analysis complete for ${meta.sample_id}"
    """
}
