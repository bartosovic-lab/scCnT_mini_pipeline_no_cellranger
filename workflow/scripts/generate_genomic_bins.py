import sys
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--chromsizes", help="Input file with chromosome lengths")
parser.add_argument("--out", help="Output file with genomic bins")
parser.add_argument("--bin", type=int, help="Size of genomic bins") 
args = parser.parse_args()


with open(args.out, "w") as f:
    for line in open(args.chromsizes):
        line = line.strip().split()
        chrom = line[0]
        length = int(line[1])
        for i in range(0,length,args.bin):
            f.write("{}\t{}\t{}\n".format(chrom,i,i+args.bin))
    