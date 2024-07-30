import pysam
import argparse
import gzip
from collections import Counter


parser = argparse.ArgumentParser(description='Summarise peak barcodes')
parser.add_argument('--bam', help='BAM file')
parser.add_argument('--out', help='Output file')
parser.add_argument('--whitelist', help='Whitelist file')
args = parser.parse_args()

def revcompl(seq):
    return seq.translate(str.maketrans('ATCG', 'TAGC'))[::-1]

# Read the whitelist
whitelist = set()
with gzip.open(args.whitelist,'rt') as f:
    for line in f:
        whitelist.add(line.strip())



with pysam.AlignmentFile(args.bam, 'rb') as bam:
    i = 0
    peak_barcodes = Counter()
    for read in bam:
        i+= 1
        peak_barcodes[read.get_tag('CB')] += 1

with open(args.out, 'w') as f:
    for key, value in peak_barcodes.items():
        if revcompl(key) in whitelist:
            f.write('{}\t{}\n'.format(key, value))
