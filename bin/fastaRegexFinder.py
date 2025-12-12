#!/usr/bin/env python3

"""
Usage:
    python fastaRegexFinder.py -f <fasta_file> [-r <regex>] [options]

Arguments:
    -f, --fasta    : Input FASTA file (use '-' for stdin)
    -r, --regex    : Regex pattern to search (default: G-quadruplex pattern)
    -m, --matchcase: Match case (default: case-insensitive)
    --noreverse    : Don't search reverse complement
    --maxstr       : Max match length to report (default: 10000)
    -s, --seqnames : List of sequence names to search
    -q, --quiet    : Suppress progress messages

Output:
    BED-format output with columns:
    chrom, start, end, match_id, length, strand, sequence

Example:
    python fastaRegexFinder.py -f reference.fa -r 'ACTG'
    python fastaRegexFinder.py -f reference.fa --seqnames chr1 chr2
"""

import re
import sys
import argparse
import operator

VERSION = '0.2.0'

parser = argparse.ArgumentParser(
    description='Search FASTA file for regex matches and output BED coordinates.',
    formatter_class=argparse.RawTextHelpFormatter
)

parser.add_argument('--fasta', '-f',
                   type=str,
                   help='Input fasta file to search. Use \'-\' for stdin.',
                   required=True)

parser.add_argument('--regex', '-r',
                   type=str,
                   help='Regex to search (default: G-quadruplex pattern)',
                   default=r'([gG]{3,}\w{1,7}){3,}[gG]{3,}')

parser.add_argument('--matchcase', '-m',
                   action='store_true',
                   help='Match case (default: case-insensitive)')

parser.add_argument('--noreverse',
                   action='store_true',
                   help='Do not search reverse complement')

parser.add_argument('--maxstr',
                   type=int,
                   required=False,
                   default=10000,
                   help='Maximum match length to report (default: 10000)')

parser.add_argument('--seqnames', '-s',
                   type=str,
                   nargs='+',
                   default=[None],
                   required=False,
                   help='List of sequence names to search')

parser.add_argument('--quiet', '-q',
                   action='store_true',
                   help='Suppress progress messages')

parser.add_argument('--version', '-v', action='version', version='%(prog)s ' + VERSION)

args = parser.parse_args()

# Set case sensitivity flag
if args.matchcase:
    flag = 0
else:
    flag = re.IGNORECASE


def sort_table(table, cols):
    """Sort a table (list of lists) by multiple columns."""
    for col in reversed(cols):
        table = sorted(table, key=operator.itemgetter(col))
    return table


def trimMatch(x, n):
    """Trim match string to max length n with notation."""
    if len(x) > n and n is not None:
        m = x[0:n] + '[' + str(n) + ',' + str(len(x)) + ']'
    else:
        m = x
    return m


def revcomp(x):
    """Reverse complement a DNA sequence. Handles ambiguity codes."""
    compdict = {
        'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
        'R': 'Y', 'Y': 'R', 'S': 'W', 'W': 'S',
        'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H',
        'H': 'D', 'V': 'B', 'N': 'N',
        'a': 't', 'c': 'g', 'g': 'c', 't': 'a',
        'r': 'y', 'y': 'r', 's': 'w', 'w': 's',
        'k': 'm', 'm': 'k', 'b': 'v', 'd': 'h',
        'h': 'd', 'v': 'b', 'n': 'n'
    }
    xrc = []
    for n in x:
        xrc.append(compdict.get(n, n))
    xrc = ''.join(xrc)[::-1]
    return xrc


def chrom_name(header):
    """Extract chromosome name from FASTA header."""
    if not header.startswith('>'):
        raise Exception('FASTA header does not start with ">":\n%s' % header)
    chr_name = re.sub(r'^>\s*', '', header)
    chr_name = re.sub(r'\s.*', '', chr_name)
    return chr_name


# Main processing
psq_re_f = re.compile(args.regex, flags=flag)

if args.fasta != '-':
    ref_seq_fh = open(args.fasta)
else:
    ref_seq_fh = sys.stdin

ref_seq = []
line = ref_seq_fh.readline()
chr_current = chrom_name(line)
line = ref_seq_fh.readline()
gquad_list = []
eof = False

while True:
    if not args.quiet:
        sys.stderr.write('Processing %s\n' % chr_current)
    
    while not line.startswith('>'):
        ref_seq.append(line.strip())
        line = ref_seq_fh.readline()
        if line == '':
            eof = True
            break
    
    ref_seq = ''.join(ref_seq)
    
    if args.seqnames == [None] or chr_current in args.seqnames:
        # Search forward strand
        for m in re.finditer(psq_re_f, ref_seq):
            matchstr = trimMatch(m.group(0), args.maxstr)
            quad_id = chr_current + '_' + str(m.start()) + '_' + str(m.end()) + '_for'
            gquad_list.append([chr_current, m.start(), m.end(), quad_id, len(m.group(0)), '+', matchstr])
        
        # Search reverse strand
        if not args.noreverse:
            ref_seq_rc = revcomp(ref_seq)
            seqlen = len(ref_seq_rc)
            for m in re.finditer(psq_re_f, ref_seq_rc):
                matchstr = trimMatch(revcomp(m.group(0)), args.maxstr)
                mstart = seqlen - m.end()
                mend = seqlen - m.start()
                quad_id = chr_current + '_' + str(mstart) + '_' + str(mend) + '_rev'
                gquad_list.append([chr_current, mstart, mend, quad_id, len(m.group(0)), '-', matchstr])
        
        # Sort and output
        gquad_sorted = sort_table(gquad_list, (1, 2, 3))
        gquad_list = []
        for xline in gquad_sorted:
            xline = '\t'.join([str(x) for x in xline])
            print(xline)
    
    if eof:
        break
    
    chr_current = chrom_name(line)
    ref_seq = []
    line = ref_seq_fh.readline()
    if line == '':
        break

sys.exit()
