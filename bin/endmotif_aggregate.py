#!/usr/bin/env python3

"""
Usage:
    python endmotif_aggregate.py <input_directory> <output_file>

Arguments:
    input_directory : Directory containing *_endmotif_ordered.csv files
    output_file     : Output CSV file path

Output:
    CSV file with motif column and percentage columns for each sample
"""

import os
import pandas as pd
import sys


def aggregate_endmotif(input_dir, output_file):
    """Aggregate end motif percentages from multiple samples into a matrix."""
    
    # Find all endmotif ordered CSV files
    endmotif_files = sorted([f for f in os.listdir(input_dir) 
                             if f.endswith('_endmotif_ordered.csv')])
    
    if not endmotif_files:
        print("No end motif files found")
        pd.DataFrame().to_csv(output_file, index=False)
        return
    
    # Read reference file for motif column
    ref_df = pd.read_csv(os.path.join(input_dir, endmotif_files[0]))
    
    # Initialize result with motif column
    if 'Motif' in ref_df.columns:
        result_df = ref_df[['Motif']].copy()
    elif ref_df.shape[1] >= 1:
        result_df = pd.DataFrame({'Motif': ref_df.iloc[:, 0]})
    else:
        result_df = pd.DataFrame()
    
    # Add percentage column from each sample
    for filename in endmotif_files:
        filepath = os.path.join(input_dir, filename)
        sample_name = filename.replace('_endmotif_ordered.csv', '')
        
        try:
            df = pd.read_csv(filepath)
            
            # Extract percentage column
            if 'Percentage of reads' in df.columns:
                pct_col = df['Percentage of reads']
            elif df.shape[1] >= 2:
                pct_col = df.iloc[:, 1]
            else:
                pct_col = df.iloc[:, -1]
            
            result_df[sample_name] = pct_col.values
            
        except Exception as e:
            print(f"Error processing {filename}: {e}")
    
    # Ensure proper CSV formatting
    result_df.to_csv(output_file, index=False, float_format='%.6f')
    
    print(f"End motif aggregate written to {output_file}")
    print(f"Samples: {len(endmotif_files)}, Motifs: {len(result_df)}")


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    
    aggregate_endmotif(sys.argv[1], sys.argv[2])


if __name__ == "__main__":
    main()
