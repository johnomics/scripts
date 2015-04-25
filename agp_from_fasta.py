#!/usr/bin/env python3

# John Davey jd626@cam.ac.uk
# Begun 29 August 2014

import argparse
import re

parser = argparse.ArgumentParser(description='''Generate AGP file from one-line FASTA file.
    (written for Hmel haplotype scaffolds))
    
    -f FASTA file''')

parser.add_argument('-f', '--fastafile', type=str, required=True)

args=parser.parse_args()

header_re = re.compile(r'^>(.+)$')
seq_re = re.compile(r'([ACGT]+)?(N+)?', re.IGNORECASE)

def make_part(seq, scaffold, pos, partnum, type):
    length = len(seq)
    end = pos + length - 1
    part = '\t'.join([scaffold, str(pos), str(end), str(partnum), type])

    if type is 'W':
        contig_name = scaffold+'.'+str(partnum)
        part = '\t'.join([part, contig_name, str(1), str(length), '+\n'])
    if type is 'N':
        part = '\t'.join([part, str(length), 'fragment', 'yes\n'])

    pos += len(seq)
    partnum += 1
    return(part, pos, partnum)

try:
    with open(args.fastafile, 'r') as fasta:
        agpname = args.fastafile + ".agp"
        with open(agpname, 'w') as agp:

            scaffold = ""
            partnum = 1
            pos = 1
            for line in fasta:
                line = line.rstrip()
                header_found = header_re.match(line)
                if header_found:
                    scaffold = header_found.group(1)
                    partnum = 1
                    pos = 1
                else:
                    for m in seq_re.finditer(line):
                        (seq, ns) = m.groups()
                        if seq:
                            (part, pos, partnum) = make_part(seq, scaffold, pos, partnum, 'W')
                            agp.write(part)
                        if ns:
                            (part, pos, partnum) = make_part(ns, scaffold, pos, partnum, 'N')
                            agp.write(part)

except IOError:
        print("Cannot open file " + args.fastafile + "!")

    