import argparse
import gzip
import os
import time
import sys

"""
Marek Bartosovic
Copyright 2024

After trimming R1 and R3, some reads are filtered out. 

This script filters R2 reads (which are not trimmed) according to R1 reads.

Resulting files should have the same reads accross R1, R2 and R3.

"""
parser = argparse.ArgumentParser(description='Filter R2 according to R1')
parser.add_argument('--fastq_R1', type=str, help='R1 file',required=True)
parser.add_argument('--fastq_R2', type=str, help='R2 file',required=True)
parser.add_argument('--out', type=str, help='Filtered R2 file',required=True)
parser.add_argument('--threads', type=int, help='Number of threads',required=False, default=8)
args = parser.parse_args()

def open_file(file):
    if file.endswith('.gz'):
        return gzip.open(file, 'rt')
    if file.endswith('.fastq') or file.endswith('.fq'):
        return open(file, 'r')
    else:
        raise ValueError('File format not recognized')

def read_fastq(f):
    name = f.readline().strip()
    seq  = f.readline().strip()
    plus = f.readline().strip()
    qual = f.readline().strip()
    return name, seq, plus, qual

start = time.time()
sys.stderr.write('Processing the fastq files\n')
i = 0

with open_file(args.fastq_R1) as f, open_file(args.fastq_R2) as g, open(args.out.replace('.gz',''), 'wt') as h:
    name1, seq1, plus1, qual1 = read_fastq(f)        
    name2, seq2, plus2, qual2 = read_fastq(g)
    while True:    
        # If names match write R2 to output file and read next R1 and R2
        i += 1
        if name1.split()[0] == name2.split()[0]:
            h.write('{}\n{}\n{}\n{}\n'.format(name2, seq2, plus2, qual2))
            name1, seq1, plus1, qual1 = read_fastq(f)
            name2, seq2, plus2, qual2 = read_fastq(g) 
        # If names don't match read next R2, which might match with the current R1, od not read new R1
        else:
            name2, seq2, plus2, qual2 = read_fastq(g)
        # If R2 is empty break
        if not name2 or not name1:
            break
        if i % 10000000 == 0:
            sys.stderr.write('Processed {} reads\n'.format(i))

sys.stderr.write('Compressing the output file using bgzip\n')
os.system('bgzip -@ {} {}'.format(args.threads, args.out.replace('.gz','')))
end=time.time()
print('Time elapsed: {}'.format(end-start))
    
