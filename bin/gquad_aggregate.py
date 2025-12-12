#!/usr/bin/env python3

"""
Usage:
    python gquad_aggregate.py <input_directory> <output_prefix>

Arguments:
    input_directory : Directory containing *_gquad_count.csv and *_gquad_total.csv files
    output_prefix   : Prefix for output files

Output:
    Three CSV files:
    - {prefix}_counts.csv  : G-quad counts per sample
    - {prefix}_totals.csv  : Total reads per sample
    - {prefix}_ratios.csv  : Ratio (count/total) per sample
"""

import os
import pandas as pd
import sys


def aggregate_gquad(input_dir, output_prefix):
    """Aggregate G-quadruplex counts and totals from multiple samples."""
    
    # Find all count and total files
    count_files = sorted([f for f in os.listdir(input_dir) 
                          if f.endswith('_gquad_count.csv')])
    total_files = sorted([f for f in os.listdir(input_dir) 
                          if f.endswith('_gquad_total.csv')])
    
    if not count_files or not total_files:
        print("No G-quad files found")
        pd.DataFrame(columns=['sample_id', 'gquad_count']).to_csv(
            f"{output_prefix}_counts.csv", index=False)
        pd.DataFrame(columns=['sample_id', 'total_reads']).to_csv(
            f"{output_prefix}_totals.csv", index=False)
        pd.DataFrame(columns=['sample_id', 'gquad_ratio']).to_csv(
            f"{output_prefix}_ratios.csv", index=False)
        return
    
    # Collect data from each sample
    counts_data = []
    totals_data = []
    ratios_data = []
    
    for count_file in count_files:
        sample_name = count_file.replace('_gquad_count.csv', '')
        total_file = f"{sample_name}_gquad_total.csv"
        
        if total_file not in total_files:
            print(f"Warning: No matching total file for {sample_name}")
            continue
        
        try:
            # Read count file
            count_df = pd.read_csv(os.path.join(input_dir, count_file))
            if 'gquad_count' in count_df.columns:
                count_val = count_df['gquad_count'].iloc[0]
            else:
                count_val = count_df.iloc[0, 1] if count_df.shape[1] > 1 else count_df.iloc[0, 0]
            
            # Read total file
            total_df = pd.read_csv(os.path.join(input_dir, total_file))
            if 'total_reads' in total_df.columns:
                total_val = total_df['total_reads'].iloc[0]
            else:
                total_val = total_df.iloc[0, 1] if total_df.shape[1] > 1 else total_df.iloc[0, 0]
            
            # Calculate ratio
            ratio_val = count_val / total_val if total_val > 0 else 0
            
            counts_data.append({'sample_id': sample_name, 'gquad_count': int(count_val)})
            totals_data.append({'sample_id': sample_name, 'total_reads': int(total_val)})
            ratios_data.append({'sample_id': sample_name, 'gquad_ratio': ratio_val})
            
        except Exception as e:
            print(f"Error processing {sample_name}: {e}")
    
    # Create DataFrames and save
    counts_df = pd.DataFrame(counts_data)
    totals_df = pd.DataFrame(totals_data)
    ratios_df = pd.DataFrame(ratios_data)
    
    # Ensure proper CSV formatting
    counts_df.to_csv(f"{output_prefix}_counts.csv", index=False)
    totals_df.to_csv(f"{output_prefix}_totals.csv", index=False)
    ratios_df.to_csv(f"{output_prefix}_ratios.csv", index=False, float_format='%.10f')
    
    print(f"G-quad aggregates written to {output_prefix}_*.csv")
    print(f"Samples: {len(counts_data)}")


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    
    aggregate_gquad(sys.argv[1], sys.argv[2])


if __name__ == "__main__":
    main()
