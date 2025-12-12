#!/usr/bin/env python3

"""
Usage:
    python insertsize_aggregate.py <input_directory> <output_file>

Arguments:
    input_directory : Directory containing *_insertsize.csv files
    output_file     : Output CSV file path

Output:
    CSV file with size column and count columns for each sample
"""

import os
import pandas as pd
import sys


def aggregate_insertsize(input_dir, output_file):
    """Aggregate insert size distributions from multiple samples into a matrix."""
    
    # Find all insertsize CSV files
    insertsize_files = sorted([f for f in os.listdir(input_dir) 
                               if f.endswith('_insertsize.csv')])
    
    if not insertsize_files:
        print("No insert size files found")
        pd.DataFrame().to_csv(output_file, index=False)
        return
    
    # Create size column (1-300 by default, or from first file)
    ref_df = pd.read_csv(os.path.join(input_dir, insertsize_files[0]))
    
    if 'size' in ref_df.columns:
        result_df = pd.DataFrame({'size': ref_df['size']})
    else:
        # Default to 1-300
        result_df = pd.DataFrame({'size': range(1, 301)})
    
    # Add count column from each sample
    for filename in insertsize_files:
        filepath = os.path.join(input_dir, filename)
        sample_name = filename.replace('_insertsize.csv', '')
        
        try:
            df = pd.read_csv(filepath)
            
            # Extract count column
            if 'count' in df.columns:
                count_col = df['count']
            elif df.shape[1] >= 2:
                count_col = df.iloc[:, 1]
            else:
                count_col = df.iloc[:, 0]
            
            # Ensure correct length (pad or truncate)
            expected_len = len(result_df)
            if len(count_col) < expected_len:
                count_col = pd.concat([count_col, 
                    pd.Series([0] * (expected_len - len(count_col)))]).reset_index(drop=True)
            elif len(count_col) > expected_len:
                count_col = count_col.iloc[:expected_len]
            
            result_df[sample_name] = count_col.values
            
        except Exception as e:
            print(f"Error processing {filename}: {e}")
    
    # Ensure proper CSV formatting (integers for counts)
    result_df.to_csv(output_file, index=False)
    
    print(f"Insert size aggregate written to {output_file}")
    print(f"Samples: {len(insertsize_files)}, Sizes: {len(result_df)}")


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    
    aggregate_insertsize(sys.argv[1], sys.argv[2])


if __name__ == "__main__":
    main()
