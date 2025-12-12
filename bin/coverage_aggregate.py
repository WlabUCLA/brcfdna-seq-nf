#!/usr/bin/env python3

"""
Usage:
    python coverage_aggregate.py <input_directory> <output_original> <output_normalized>

Arguments:
    input_directory   : Directory containing *_multibam.csv files
    output_original   : Output CSV with raw coverage values
    output_normalized : Output CSV with per-sample normalized values

Output:
    Two CSV files with coverage data aggregated across samples
"""

import os
import pandas as pd
import numpy as np
import sys


def aggregate_coverage(input_dir, output_original, output_normalized):
    """Aggregate coverage data from multiple samples into matrices."""
    
    # Find all multibam CSV files
    coverage_files = sorted([f for f in os.listdir(input_dir) 
                             if f.endswith('_multibam.csv')])
    
    if not coverage_files:
        print("No coverage files found")
        pd.DataFrame().to_csv(output_original, index=False)
        pd.DataFrame().to_csv(output_normalized, index=False)
        return
    
    # Read reference file for bin information
    ref_df = pd.read_csv(os.path.join(input_dir, coverage_files[0]))
    
    # Initialize result with bin columns (chr, start, end if present)
    bin_cols = []
    for col in ['chr', 'chrom', 'chromosome', 'start', 'end', 'bin']:
        if col in ref_df.columns:
            bin_cols.append(col)
    
    if bin_cols:
        result_df = ref_df[bin_cols].copy()
    else:
        result_df = pd.DataFrame({'bin': range(len(ref_df))})
    
    # Collect coverage data from each sample
    sample_data = {}
    
    for filename in coverage_files:
        filepath = os.path.join(input_dir, filename)
        sample_name = filename.replace('_multibam.csv', '')
        
        try:
            df = pd.read_csv(filepath)
            
            # Find coverage column (usually last numeric column)
            numeric_cols = df.select_dtypes(include=[np.number]).columns
            if len(numeric_cols) > 0:
                coverage_col = df[numeric_cols[-1]]
            else:
                coverage_col = df.iloc[:, -1]
            
            sample_data[sample_name] = coverage_col.values
            
        except Exception as e:
            print(f"Error processing {filename}: {e}")
    
    # Add sample columns to result
    for sample_name, coverage in sample_data.items():
        expected_len = len(result_df)
        if len(coverage) < expected_len:
            coverage = np.concatenate([coverage, np.zeros(expected_len - len(coverage))])
        elif len(coverage) > expected_len:
            coverage = coverage[:expected_len]
        
        result_df[sample_name] = coverage
    
    # Save original values
    result_df.to_csv(output_original, index=False, float_format='%.6f')
    
    # Create normalized version (per-sample normalization)
    normalized_df = result_df.copy()
    sample_cols = [c for c in normalized_df.columns if c not in bin_cols]
    
    for col in sample_cols:
        col_sum = normalized_df[col].sum()
        if col_sum > 0:
            normalized_df[col] = normalized_df[col] / col_sum
    
    normalized_df.to_csv(output_normalized, index=False, float_format='%.10f')
    
    print(f"Coverage aggregates written to:")
    print(f"  Original: {output_original}")
    print(f"  Normalized: {output_normalized}")
    print(f"Samples: {len(sample_data)}, Bins: {len(result_df)}")


def main():
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)
    
    aggregate_coverage(sys.argv[1], sys.argv[2], sys.argv[3])


if __name__ == "__main__":
    main()
