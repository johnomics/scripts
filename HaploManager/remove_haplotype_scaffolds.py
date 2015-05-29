#!/usr/bin/python3 -u

import argparse
from Bio import SeqIO

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='''Output TSV for new genome
        -f fasta
        -b broken_haplotypes
        -p prefix
        -o output
        ''')

    parser.add_argument('-f', '--fasta', type=str, required=True)
    parser.add_argument('-b', '--broken', type=str, required=False)
    parser.add_argument('-p', '--prefix', type=str, required=True)
    parser.add_argument('-o', '--output', type=str, required=True)

    args = parser.parse_args()
    
    broken = {}
    if args.broken:
        brokefile = open(args.broken, 'r')
        for line in brokefile:
            scaffold = line.rstrip()
            broken[scaffold] = 1
    
    fin = open(args.fasta, 'r')
    fout = open(args.output, 'w')
    for seq in SeqIO.parse(fin, "fasta"):
        if seq.name in broken or args.prefix in seq.name:
            continue
        SeqIO.write(seq, fout, "fasta")
        