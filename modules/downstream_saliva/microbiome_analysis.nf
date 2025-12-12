/*
================================================================================
    MICROBIOME ANALYSIS MODULE (MetaPhlAn 4)
================================================================================
    Performs microbiome profiling on unmapped reads using MetaPhlAn 4.
    Identifies bacterial, archaeal, viral, and eukaryotic species in saliva samples.
    
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
    tuple val(meta), path(bam)
    
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
    echo "Database: ${params.metaphlan_db ?: 'default (auto-detect)'}"
    
    # Extract unmapped reads from BAM (flag 4 = unmapped)
    echo "Extracting unmapped reads..."
    samtools view -f 4 -b ${bam} | \\
        samtools fastq -@ ${task.cpus} - | \\
        gzip > ${meta.sample_id}_unmapped.fastq.gz
    
    # Count unmapped reads
    UNMAPPED_COUNT=\$(zcat ${meta.sample_id}_unmapped.fastq.gz | wc -l)
    UNMAPPED_COUNT=\$((UNMAPPED_COUNT / 4))
    echo "Unmapped reads extracted: \${UNMAPPED_COUNT}"
    
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
        echo "=== Summary ==="
        head -20 ${meta.sample_id}_mpa.csv
        
    else
        echo "WARNING: No unmapped reads found in ${bam}"
        echo "clade_name,relative_abundance" > ${meta.sample_id}_mpa.csv
        echo "UNKNOWN,100.0" >> ${meta.sample_id}_mpa.csv
    fi
    
    echo "Microbiome analysis complete for ${meta.sample_id}"
    """
}

/*
 * Optional: Database installation process
 * Run this once to set up the MetaPhlAn database
 */
process METAPHLAN_INSTALL_DB {
    label 'process_high'
    storeDir "${params.metaphlan_db ?: 'metaphlan_db'}"
    
    output:
    path "mpa_vJan21_CHOCOPhlAnSGB_202103/*", emit: db_files
    
    script:
    def db_path = params.metaphlan_db ?: 'metaphlan_db'
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    echo "Installing MetaPhlAn database to ${db_path}..."
    echo "This may take 10-30 minutes depending on network speed."
    
    metaphlan --install --bowtie2db ${db_path}
    
    echo "Database installation complete."
    echo "Files installed:"
    ls -lh ${db_path}/
    """
}
