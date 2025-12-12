#!/usr/bin/env python3

"""
Usage:
    python fragmentomic_ratio.py <num_chromosomes> <window_length> <num_bins> <min_short> <max_short> <max_long> <input_sam>

Arguments:
    num_chromosomes  : Number of chromosomes (typically 25)
    window_length    : Bin window size in bp (typically 1000000)
    num_bins         : Number of bins per chromosome (typically 300)
    min_short        : Minimum fragment length for short category (e.g., 25)
    max_short        : Maximum fragment length for short category (e.g., 100)
    max_long         : Maximum fragment length for long category (e.g., 250)
    input_sam        : Input SAM file

Output:
    Creates <input_sam>_fragmentomic_ratio.csv with short/long fragment ratios per genomic bin
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for server use
import matplotlib.pyplot as plt
import sys


class frag():
    def __init__(self, numberofchr, windows_length, numberofbin, minsizeA, maxsizeA, maxsizeB, inputfile, outputfile):
       
        self.numberofchr = numberofchr
        self.windows_length = windows_length
        self.numberofbin = numberofbin
        self.inputfile = inputfile
        self.outputfile = outputfile
        self.minsizeA = minsizeA
        self.maxsizeA = maxsizeA
        self.maxsizeB = maxsizeB
        self.allinfoA = np.zeros((self.numberofchr, self.numberofbin))
        self.allinfoB = np.zeros((self.numberofchr, self.numberofbin))
        
        self.readfile(self.inputfile)
        
        print(self.allinfoA)
        print(self.allinfoB)
        self.writeinfile(self.outputfile)

    def readfile(self, filename):
        with open(filename) as fp:
            line = fp.readline()
            cnt = 1
            while line:
                line = str(line).strip().split()
                cnt += 1   
                self.analyse_line(line, self.minsizeA, self.maxsizeA, self.maxsizeB)
                line = fp.readline()
                
    def analyse_lines(self, list_line):
        pass
    
    def analyse_line(self, line, minsizeA, maxsizeA, maxsizeB):
        if len(line) > 10:
            chrom = str(line[2].replace('chr', ''))
           
            if 'X' in chrom or 'x' in chrom:
                chrom = 23
            elif 'Y' in chrom or 'y' in chrom:
                chrom = 24
            else:
                try:
                    chrom = int(chrom)
                except ValueError:
                    return  # Skip non-standard chromosomes
            
            if chrom >= self.numberofchr:
                return  # Skip if chromosome index out of range
            
            bin_index = int(int(line[3]) / self.windows_length)
            
            if bin_index >= self.numberofbin:
                return  # Skip if bin index out of range
            
            lengthof = len(line[9])

            if lengthof >= minsizeA and lengthof < maxsizeA:
                self.allinfoA[chrom, bin_index] += 1
                
            elif lengthof >= maxsizeA and lengthof < maxsizeB:
                self.allinfoB[chrom, bin_index] += 1
           
    def writeinfile(self, outputfile):
        with open(outputfile, 'w') as f:
            f.write("chromosome,bin_index,reads_short,reads_long,normalized_short,normalized_long,ratio\n")
            for x in range(1, self.numberofchr):  # Start from chromosome 1
                for y in range(self.numberofbin):
                    num_reads_A = self.allinfoA[x, y]
                    num_reads_B = self.allinfoB[x, y]
                    total = num_reads_A + num_reads_B
                    
                    norm_A = num_reads_A / total if total != 0 else 0
                    norm_B = num_reads_B / total if total != 0 else 0
                    ratio = norm_A / norm_B if norm_B != 0 else 0
                    
                    line = f"{x},{y},{int(num_reads_A)},{int(num_reads_B)},{norm_A:.6f},{norm_B:.6f},{ratio:.6f}\n"
                    f.write(line)


def main(argv):
    if len(argv) < 7:
        print(__doc__)
        sys.exit(1)
        
    numberofchr = int(argv[0])
    windows_length = int(argv[1])
    numberofbin = int(argv[2])
    minsizeA = int(argv[3])
    maxsizeA = int(argv[4])
    maxsizeB = int(argv[5])
    inputfile = argv[6]
    outputfile = inputfile + '_fragmentomic_ratio.csv'
    
    f = frag(numberofchr, windows_length, numberofbin, minsizeA, maxsizeA, maxsizeB, inputfile, outputfile)


if __name__ == "__main__":
    main(sys.argv[1:])
