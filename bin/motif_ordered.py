#!/usr/bin/env python3

"""
Usage:
    python motif_ordered.py <input_csv> <output_csv>

Arguments:
    input_csv   : Input CSV with Motif, Frequency, Percentage columns
    output_csv  : Output CSV with all 256 4-mer motifs (missing filled with 0)

Output:
    CSV file with all ACGT 4-mer combinations ordered alphabetically
"""

import pandas as pd
import itertools
import sys


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    
    input_csv = sys.argv[1]
    output_csv = sys.argv[2]
    
    # Read the input CSV file
    data = pd.read_csv(input_csv)
    
    # Generate all possible combinations of ACTG
    combinations = [''.join(comb) for comb in itertools.product('acgt', repeat=4)]
    
    # Sort the data by the "Motif" column
    data = data.sort_values(by='Motif')
    
    # Build result as a list of dictionaries (pandas 2.x compatible)
    result_rows = []
    
    # Fill in missing combinations with 0s
    for combo in combinations:
        if combo in data['Motif'].values:
            row = data[data['Motif'] == combo].iloc[0]
            result_rows.append({
                'Motif': row['Motif'],
                'Frequency': row['Frequency'],
                'Percentage of reads': row['Percentage of reads']
            })
        else:
            result_rows.append({
                'Motif': combo,
                'Frequency': 0,
                'Percentage of reads': 0
            })
    
    # Create DataFrame from list of dictionaries
    result = pd.DataFrame(result_rows)
    
    # Write the result to a new CSV file
    result.to_csv(output_csv, index=False)
    
    print(f"Output has been written to {output_csv}")


if __name__ == "__main__":
    main()
