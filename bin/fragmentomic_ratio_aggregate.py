#!/usr/bin/env python3

"""
Usage:
    python fragmentomic_ratio_aggregate.py <input_directory> <output_file>

Arguments:
    input_directory : Directory containing *_fragmentomic_ratio.csv files
    output_file     : Output CSV file path

Output:
    CSV file with chromosome, bin_index, and ratio columns for each sample
"""

import os
import pandas as pd
import sys


def aggregate_fragmentomic_ratio(input_dir, output_file):
    """Aggregate fragmentomic ratio values from multiple samples into a matrix."""
    
    # Find all fragmentomic ratio CSV files
    frag_files = sorted([f for f in os.listdir(input_dir) if f.endswith('_fragmentomic_ratio.csv')])
    
    if not frag_files:
        print("No fragmentomic ratio files found")
        pd.DataFrame().to_csv(output_file, index=False)
        return
    
    # Read reference file for chromosome and bin_index columns
    ref_df = pd.read_csv(os.path.join(input_dir, frag_files[0]))
    
    # Initialize result with reference columns
    if 'chromosome' in ref_df.columns and 'bin_index' in ref_df.columns:
        result_df = ref_df[['chromosome', 'bin_index']].copy()
    else:
        result_df = pd.DataFrame()
    
    # Add ratio column from each sample
    for filename in frag_files:
        filepath = os.path.join(input_dir, filename)
        sample_name = filename.replace('_fragmentomic_ratio.csv', '')
        
        try:
            df = pd.read_csv(filepath)
            
            # Extract ratio column
            if 'ratio' in df.columns:
                ratio_col = df['ratio']
            elif df.shape[1] >= 3:
                ratio_col = df.iloc[:, 2]
            else:
                ratio_col = df.iloc[:, -1]
            
            result_df[sample_name] = ratio_col.values
            
        except Exception as e:
            print(f"Error processing {filename}: {e}")
    
    # Ensure proper CSV formatting
    result_df.to_csv(output_file, index=False, float_format='%.6f')
    
    print(f"Fragmentomic ratio aggregate written to {output_file}")
    print(f"Samples: {len(frag_files)}, Rows: {len(result_df)}")


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    
    aggregate_fragmentomic_ratio(sys.argv[1], sys.argv[2])


if __name__ == "__main__":
    main()
