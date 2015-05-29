#!/usr/bin/python3 -u

import argparse
from Bio import SeqIO

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='''Output TSV for new genome
        -f fasta
        -o output
        ''')

    parser.add_argument('-f', '--fasta', type=str, required=True)
    parser.add_argument('-o', '--output', type=str, required=True)

    args = parser.parse_args()
    
    tsv = open(args.output, 'w')
    f = open(args.fasta, 'r')
    for seq in SeqIO.parse(f, "fasta"):
        end = len(seq)
        tsv.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(seq.name, 1, end, seq.name, 1, end, 1, "active"))