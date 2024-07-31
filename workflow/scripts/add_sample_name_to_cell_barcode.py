import argparse
import os
import gzip


parser = argparse.ArgumentParser()
parser.add_argument('--fragments', type=str, required=True)
parser.add_argument('--output', type=str, required=True)
parser.add_argument('--sample_name', type=str, required=True)
parser.add_argument('--threads', type=int, default=8)
args = parser.parse_args()


def open_file(file_path):
    if file_path.endswith('.gz'):
        return gzip.open(file_path, 'rt')
    if file_path.endswith('.tsv') or file_path.endswith('.txt') or file_path.endswith('.bed'):
        return open(file_path, 'r')
    return None

def main():
    with open_file(args.fragments) as f:
        with open(args.output.replace('.gz', ''), 'w') as o:
            for line in f:
                line = line.strip().split('\t')
                o.write('{}\t{}\t{}\t{}_{}\t{}\n'.format(line[0], line[1], line[2], line[3], args.sample_name, line[4]))
    os.system('bgzip -@ {} -f {}'.format(args.threads, args.output.replace('.gz', '')))


if __name__ == '__main__':
    main()