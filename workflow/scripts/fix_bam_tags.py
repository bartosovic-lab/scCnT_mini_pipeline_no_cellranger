import pysam
import argparse
import sys 

parser = argparse.ArgumentParser(description='Fix BAM tags')
parser.add_argument('--bam', type=str, help='BAM file',required=True)
parser.add_argument('--out', type=str, help='Output BAM file',required=True)
parser.add_argument('--sep', type=str, help='Separator',required=False, default='|')
args = parser.parse_args()

def remove_CR_CY(string):
    string = string.replace('CR:Z:', '')
    string = string.replace('CY:Z:', '')
    return string


def main():
    with pysam.AlignmentFile(args.bam, 'rb') as f, pysam.AlignmentFile(args.out, 'wb', template=f) as g:
        for read in f:
            name = read.query_name.split(args.sep)
            CR_tag = remove_CR_CY(name[0])
            CY_tag = remove_CR_CY(name[1])
            if len(CR_tag) != 16 or len(CY_tag) != 16:
                sys.stderr.write('Barcode length is not 16 for read: {}\n'.format(name))
                continue
            read.set_tag('CR', CR_tag)
            read.set_tag('CY', CY_tag)
            read.set_tag('CB', CR_tag)
            read.query_name = name[2]
            g.write(read)

if __name__ == '__main__':
    main()
