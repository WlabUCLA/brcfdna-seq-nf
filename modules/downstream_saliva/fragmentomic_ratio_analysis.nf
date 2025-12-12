/*
================================================================================
    FRAGMENTOMIC RATIO ANALYSIS MODULE
================================================================================
    Analyzes fragment length distributions in genomic bins.
    Calculates short/long fragment ratios for cancer detection signatures.
    
    Outputs:
    - Raw fragmentomic ratio results with all columns
    - Parsed results (chromosome, bin, ratio)
*/

process FRAGMENTOMIC_RATIO_ANALYSIS {
    tag "${meta.sample_id}"
    label 'process_medium'
    
    publishDir "${params.outdir}/downstream/fragmentomic_ratio/${meta.sample_id}", mode: 'copy'
    
    input:
    tuple val(meta), path(bam)
    
    output:
    tuple val(meta), path("${meta.sample_id}_fragmentomic_ratio_raw.csv"), emit: fragmentomic_raw
    tuple val(meta), path("${meta.sample_id}_fragmentomic_ratio.csv"), emit: fragmentomic_csv
    
    script:
    def num_chrom = 25
    def num_bins = 300
    def min_short = params.fragmentomic_min_short ?: 25
    def max_short = params.fragmentomic_max_short ?: 100
    def max_long = params.fragmentomic_max_long ?: 250
    def bin_size = params.fragmentomic_bin_size ?: 1000000
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    # Convert BAM to SAM for fragmentomic ratio processing
    samtools view ${bam} -o ${meta.sample_id}.sam -O sam
    
    # Run fragmentomic ratio analysis
    python3 ${projectDir}/bin/fragmentomic_ratio.py \\
        ${num_chrom} \\
        ${bin_size} \\
        ${num_bins} \\
        ${min_short} \\
        ${max_short} \\
        ${max_long} \\
        ${meta.sample_id}.sam
    
    # Rename raw output
    mv ${meta.sample_id}.sam_delfi_frag.csv ${meta.sample_id}_fragmentomic_ratio_raw.csv
    
    # Parse output - extract chromosome, bin, and ratio columns
    python3 << 'PYTHON_SCRIPT'
import csv

input_file = "${meta.sample_id}_fragmentomic_ratio_raw.csv"
output_file = "${meta.sample_id}_fragmentomic_ratio.csv"

# Columns to extract: B (1), D (3), N (13) - 0-indexed
# These are: chromosome, bin_index, ratio
columns_to_extract = [1, 3, 13]

with open(input_file, 'r', newline='') as infile, open(output_file, 'w', newline='') as outfile:
    reader = csv.reader(infile)
    writer = csv.writer(outfile)
    
    # Write header
    writer.writerow(['chromosome', 'bin_index', 'ratio'])
    
    for row in reader:
        if len(row) > max(columns_to_extract):
            extracted = [row[i] for i in columns_to_extract]
            writer.writerow(extracted)

print(f"Parsed fragmentomic ratio output written to {output_file}")
PYTHON_SCRIPT
    
    # Clean up SAM file
    rm -f ${meta.sample_id}.sam
    
    echo "Fragmentomic ratio analysis complete for ${meta.sample_id}"
    """
}
