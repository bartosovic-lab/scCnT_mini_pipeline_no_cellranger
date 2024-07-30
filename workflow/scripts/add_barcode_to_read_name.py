import argparse
import gzip
import time
import os,sys

parser = argparse.ArgumentParser(description='Add R2 barcode to read names of R1/R3')
parser.add_argument('--fastq_R1', type=str, help='R1 file',required=True)
parser.add_argument('--fastq_R2', type=str, help='R2 file',required=True)
parser.add_argument('--fastq_R3', type=str, help='R3 file',required=True)
parser.add_argument('--out_R1', type=str, help='R1 output file',required=True)
parser.add_argument('--out_R3', type=str, help='R3 output file',required=True)
parser.add_argument('--threads', type=int, help='Number of threads',required=False,default=8)
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

sys.stderr.write('Reading files\n')
start = time.time()
i=0 

with open_file(args.fastq_R1) as f1, open_file(args.fastq_R2) as f2, open_file(args.fastq_R3) as f3, open(args.out_R1.replace('.gz',''), 'wt') as g1, open(args.out_R3.replace('.gz',''), 'wt') as g2:
    name1, seq1, plus1, qual1 = read_fastq(f1)        
    name2, seq2, plus2, qual2 = read_fastq(f2)
    name3, seq3, plus3, qual3 = read_fastq(f3)
    while True:
        # If names match write R2 to output file and read next R1 and R2
        i+=1
        if name1.split()[0] == name2.split()[0] and name3.split()[0] == name2.split()[0]:
            g1.write('@CR:Z:{}|CY:Z:{}|{}\n{}\n{}\n{}\n'.format(seq2,qual2,name1.lstrip("@"), seq1, plus1, qual1))
            g2.write('@CR:Z:{}|CY:Z:{}|{}\n{}\n{}\n{}\n'.format(seq2,qual2,name3.lstrip("@"), seq3, plus3, qual3))
            name1, seq1, plus1, qual1 = read_fastq(f1)
            name2, seq2, plus2, qual2 = read_fastq(f2)
            name3, seq3, plus3, qual3 = read_fastq(f3)
        else:
            raise ValueError('Names of reads do not match')
        # If R2 is empty break
        if not name2:
            break
        if i%1000000==0:
            sys.stderr.write('Processed {} reads\n'.format(i))

sys.stderr.write('Compressing files\n')
os.system('bgzip -@ {} {} '.format(args.threads, args.out_R1.replace('.gz','')))
os.system('bgzip -@ {} {} '.format(args.threads, args.out_R3.replace('.gz','')))
end = time.time()
print("Time taken: ", end-start)