#!/usr/bin/python3 -u

import os
import sys
import pathlib
import argparse
import glob
import sqlite3 as sql
from pprint import pprint
from Bio import SeqIO
from collections import defaultdict, namedtuple
from operator import itemgetter

Haplotype = namedtuple("Haplotype", "hapname hapstart hapend hapstrand name start end strand")

class Part():
    def __init__(self, oldname, oldstart, oldend, newname, newstart, newend, strand, parttype, hap=None):
        self.oldname = oldname
        self.oldstart = int(oldstart)
        self.oldend = int(oldend)
        self.newname = newname
        self.newstart = int(newstart)
        self.newend = int(newend)
        self.strand = int(strand)
        self.parttype = parttype
        self.hap = None
        if hap:
            self.hap = Haplotype(hap[0], int(hap[1]), int(hap[2]), int(hap[3]), None, None, None, None)

    @property
    def length(self):
        return self.oldend - self.oldstart + 1
    
    def __repr__(self):
        out = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(self.oldname, self.oldstart, self.oldend, self.newname,
                                                               self.newstart, self.newend, self.strand, self.parttype)
        if self.hap:
            out += '\t{}\t{}\t{}\t{}'.format(self.hap.hapname, self.hap.hapstart, self.hap.hapend, self.hap.hapstrand)
        return out

def load_genome(mergedgenome):
    genome = defaultdict(list)
    haps = defaultdict(lambda: defaultdict(list))
    try:
        with open(mergedgenome, 'r') as g:
            for line in g:
                oldname, oldstart, oldend, newname, newstart, newend, strand, typename, *args = line.rstrip().split('\t')
                genome[newname].append(Part(oldname, oldstart, oldend, newname, newstart, newend, strand, typename, args))
                if args:
                    hapname, hapstart, hapend, hapstrand = args
                    haps[hapname][oldname].append(Haplotype(hapname, int(hapstart), int(hapend), int(hapstrand), oldname, int(oldstart), int(oldend), int(strand)))
        return genome, haps
    except IOError:
        print("Failed to load genome file {}".format(mergedgenome))
        sys.exit()
        
def get_orig_pos(pos, part):
    origpos = None
    if pos and part.newstart <= pos <= part.newend:
        if part.strand == 1:
            origpos = part.oldstart + pos - part.newstart
        elif part.strand == -1:
            origpos = part.oldstart + part.newend - pos
    return origpos

def order(start, end):
    if end is not None and start is not None and start > end:
        start, end = end, start
    return start, end

def check_haplotypes(name, origstart, origend, genome):
    if name in haps:
        for hapname in haps[name]:
            for hap in haps[name][hapname]:
                if origstart and hap.hapstart <= origstart <= hap.hapend or origend and hap.hapstart <= origend <= hap.hapend:
                    start, end = order(hap.start, hap.end)
                    print("{} {}:{} {} is haplotype of {} {}:{} {}".format(hap.name, hap.start, hap.end, hap.strand, hap.hapname, hap.hapstart, hap.hapend, hap.hapstrand))
                    print('M\t{}\t{}\t{}\tR\t{}-{}'.format(hap.name, start, end, start, end))

def load_misassemblies(misfile, genome, haps):
    try:
        with open(misfile, 'r') as m:
            for line in m:
                p, name, mistype, details = line.rstrip().split('\t')
                print(name, mistype, details)

                if name in genome:
                    start, end = None, None
                    if mistype == 'R':
                        start, end = [int(x) for x in details.split('-')]
                    elif mistype == 'B':
                        start = int(details)
                    for part in sorted(genome[name],key=lambda x:x.newstart):
                        origstart = get_orig_pos(start, part)
                        origend   = get_orig_pos(end,   part)
                        if origstart or origend:
                            print(part)
                            if mistype == 'R':
                                printstart, printend = order(origstart, origend)
                                print('M\t{}\t{}\t{}\tR\t{}-{}'.format(part.oldname, part.oldstart, part.oldend, printstart, printend))
                            elif mistype == 'B':
                                print('M\t{}\t{}\t{}\tB\t{}'.format(part.oldname, part.oldstart, part.oldend, origstart))
                            check_haplotypes(part.oldname, origstart, origend, haps)
    except IOError:
        print("Failed to load misassemblies file", misfile)
        sys.exit()

def get_args():
    parser = argparse.ArgumentParser(description='''Output misassembly positions from merged genome on original genome
    
        -t tsv
        -m misassemblies
        ''')

    parser.add_argument('-t', '--tsv', type=str, required=True)
    parser.add_argument('-m', '--misassemblies', type=str, required=True)
    return parser.parse_args()

if __name__ == '__main__':
    
    args = get_args()

    genome, haps = load_genome(args.tsv)
    
    misassemblies = load_misassemblies(args.misassemblies, genome, haps)