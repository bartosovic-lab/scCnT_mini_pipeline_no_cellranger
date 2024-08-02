import argparse
import gzip
import time
from datetime import timedelta
import os

"""
This script is useful in case there is single-modality scCUT&Tag data that was sequenced using nano-CT sequencing setup [36-8-48-36] on the same lane/flowcell

Sequencing structure of R2 in nano-CT sequencing is as follows:

PN-2000210 (10x Gel bead oligo) 5-AATGATACGGCGACCACCGAGATCTACAC-NNNNNNNNNNNNNNNN-TCGTCGGCAGCGTC-3
                                                                               5' TCGTCGGCAGCGTCTCCACGC NNNNNNNN GCGATCGAGGACGGCAGATGTGTATAAGAGACAG  3' 
                                                                single-cell barcode                 modality barcode

Sequencing structure of R2 single-cell CUT&Tag (single modality)

PN-2000210 (10x Gel bead oligo) 5'-AATGATACGGCGACCACCGAGATCTACAC-NNNNNNNNNNNNNNNN-TCGTCGGCAGCGTC-3'
                                                                               5' TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG 3' 

Therefore first 16bp are the single-cell index in case of scCUT&Tag data sequenced using nano-CT setup.
                                                                               
The script will extract the single-cell barcode from the nano-CT sequencing data and write it to a new fastq file.
"""
parser = argparse.ArgumentParser(description='Get scCT barcode from nano-CT sequencing')
parser.add_argument('--input', type=str, help='Input R2 fastq file from nano-CT sequencing [36-8-48-36]')
parser.add_argument('--output', type=str, help='Output R2 file with 16bp single-cell barcode [16bp index]')
parser.add_argument('--threads', type=int, default=8, help='Number of threads to use for bgzip compression')
args = parser.parse_args()

def open_file(file):
    if file.endswith('.gz'):
        return gzip.open(file, 'rt')
    else:
        return open(file, 'r')
    

def main():
    t0 = time.time()

    n=0
    k=100000
    barcodes = [] 
    print('*** Doing a routine check of barcode complexity')
    with open_file(args.input) as f, open(args.output.replace('.gz',''), 'w') as out:
        for line in f:
            n+=1
            if line.startswith('@'):
                
                name = line
                seq = next(f)
                plus = next(f)
                qual = next(f)
                
                out.write(name)
                out.write(seq[:16] + '\n')
                out.write(plus)
                out.write(qual[:16] + '\n')

                if n < k:
                    barcodes.append(seq[:16])
                if n == k:
                    print('{} barcodes found in the first {k} reads'.format(len(set(barcodes)), k=k))
                    if len(set(barcodes)) < 500:
                        print('*** WARNING: Low barcode diversity detected. Check the input file')
                    if len(set(barcodes)) >= 500:
                        print('*** INFO: High barcode diversity detected.')
                    del(barcodes)
                    
            # Quick check for barcode complexity
    
    os.system('bgzip -@ {} {}'.format(args.threads, args.output.replace('.gz','')))
    t1 = time.time()
    print('*** INFO: Finished')
    print('*** INFO: Time taken: {}'.format(timedelta(seconds=t1-t0)))
            
                
                
            

if __name__ == '__main__':
    main()
