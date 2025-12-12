#!/usr/bin/env python3

"""
Usage:
    python motif_count.py <input_bam> <temp_motifs_file> <output_csv>

Arguments:
    input_bam         : Input BAM file
    temp_motifs_file  : Temporary file for extracted 4-mer motifs
    output_csv        : Output CSV with motif frequencies

Output:
    CSV file with columns: Motif, Frequency, Percentage of reads
"""

import csv
import pysam
import sys


def extract_first_four_base_pairs(bam_file_path, output_file_path):
    """Extract the first 4 base pairs from every read in a BAM file."""
    with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
        with open(output_file_path, "w") as output_file:
            for read in bam_file:
                if read.query_sequence and len(read.query_sequence) >= 4:
                    first_four_base_pairs = read.query_sequence[:4]
                    output_file.write(f"{first_four_base_pairs}\n")


def count_motifs(file_path):
    """Count frequency of each 4-mer motif."""
    motif_frequency = {}
    
    try:
        with open(file_path, 'r', encoding='utf-8') as file:
            text = file.read()
            motifs = text.split()
            
            for motif in motifs:
                motif = motif.lower()
                motif_frequency[motif] = motif_frequency.get(motif, 0) + 1
        
        return motif_frequency
    
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        return None


def remove_motifs_with_n(d):
    """Remove motifs containing 'N' bases."""
    motifs_with_n = [key for key in d.keys() if 'n' in key.lower()]
    for key in motifs_with_n:
        d.pop(key)
    return d


def write_dict_to_csv_with_ratio(d, file_path):
    """Write motif frequencies to CSV with percentage column."""
    total_sum = sum(d.values())
    
    rows = [["Motif", "Frequency", "Percentage of reads"]]
    for key, value in d.items():
        value_ratio = value / total_sum if total_sum > 0 else 0
        rows.append([key, value, value_ratio])
    
    try:
        with open(file_path, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerows(rows)
        print(f"Data written to {file_path} successfully.")
    except Exception as e:
        print(f"An error occurred while writing to the CSV file: {str(e)}")


def main():
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)
    
    bam_file = sys.argv[1]
    temp_file = sys.argv[2]
    output_csv = sys.argv[3]
    
    # Extract first 4 base pairs from each read
    extract_first_four_base_pairs(bam_file, temp_file)
    
    # Count motif frequencies
    motif_frequency = count_motifs(temp_file)
    
    if motif_frequency is None:
        sys.exit(1)
    
    # Remove motifs with N's
    motifs_without_n = remove_motifs_with_n(motif_frequency)
    
    # Write to CSV
    write_dict_to_csv_with_ratio(motifs_without_n, output_csv)


if __name__ == "__main__":
    main()
