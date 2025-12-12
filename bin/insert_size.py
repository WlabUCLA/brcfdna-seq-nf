#!/usr/bin/env python3

"""
Usage:
    python insert_size.py <input_csv> <output_sorted> <output_filled>

Arguments:
    input_csv      : Input CSV file with count,size columns
    output_sorted  : Output CSV with sorted size,count columns
    output_filled  : Output CSV with all sizes 1-300 filled (missing = 0)

Output:
    Two CSV files with insert size distributions
"""

import sys


def main():
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)
    
    input_filename = sys.argv[1]
    output_filename = sys.argv[2]
    new_output = sys.argv[3]
    
    # Create a list to store the data
    data = []
    
    # Read the input file and populate the list (CSV format: count,size)
    try:
        with open(input_filename, 'r') as file:
            for line in file:
                line = line.strip()
                if not line or line.startswith('count'):  # Skip header
                    continue
                columns = line.split(",")
                if len(columns) == 2:
                    data.append((int(columns[1]), int(columns[0])))  # (size, count)
    except FileNotFoundError:
        print(f"Error: File '{input_filename}' not found.")
        sys.exit(1)
    
    # Sort the data by size
    sorted_data = sorted(data, key=lambda x: x[0])
    
    # Write the sorted data to CSV file
    with open(output_filename, 'w') as file:
        file.write("size,count\n")
        for item in sorted_data:
            file.write(f"{item[0]},{item[1]}\n")
    
    print(f"Output written to {output_filename}")
    
    # Create a dictionary to store the existing data
    data_dict = {}
    
    # Read the sorted file and populate the dictionary
    try:
        with open(output_filename, 'r') as file:
            next(file)  # Skip header
            for line in file:
                columns = line.strip().split(',')
                if len(columns) == 2:
                    key, value = columns
                    data_dict[int(key)] = int(value)
    except FileNotFoundError:
        print(f"Error: File '{output_filename}' not found.")
        sys.exit(1)
    
    # Fill in missing values with 0 from 1 to 300
    for key in range(1, 301):
        if key not in data_dict:
            data_dict[key] = 0
    
    # Sort the data by the keys in ascending order
    sorted_data = sorted(data_dict.items())
    
    # Write the sorted data to the output CSV file
    with open(new_output, 'w') as file:
        file.write("size,count\n")
        for key, value in sorted_data:
            file.write(f"{key},{value}\n")
    
    print(f"Output written to {new_output}")


if __name__ == "__main__":
    main()
