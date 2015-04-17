#!/usr/bin/python3 -u

import os
import sys
import pathlib
import argparse
import sqlite3 as sql
from collections import defaultdict
from operator import itemgetter

class Part():
    def __init__(self, oldname, oldstart, oldend, newname, newstart, newend, strand, parttype):
        self.oldname = oldname
        self.oldstart = int(oldstart)
        self.oldend = int(oldend)
        self.newname = newname
        self.newstart = int(newstart)
        self.newend = int(newend)
        self.strand = int(strand)
        self.parttype = parttype

    @property
    def length(self):
        return self.oldend - self.oldstart + 1
    
    def __repr__(self):
        return '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(self.oldname, self.oldstart, self.oldend, self.newname, self.newstart, self.newend, self.strand, self.parttype)


def load_genome(mergedgenome):
    genome = defaultdict(list)
    try:
        with open(mergedgenome, 'r') as g:
            for line in g:
                oldname, oldstart, oldend, newname, newstart, newend, strand, parttype = line.rstrip().split('\t')
                if strand == '0':
                    strand = '1'
                genome[oldname].append(Part(oldname, oldstart, oldend, newname, newstart, newend, strand, typename))
    except IOError:
        print("Failed to load genome file {}".format(mergedgenome))
        sys.exit()
    
    return genome

def get_parts(scaffold, start, end, genome):

    parts = []
    for part in genome[scaffold]:
        if part.oldstart > end or part.oldend < start:
            continue
        slice_start = max(part.oldstart, start)
        slice_end = min(part.oldend, end)
        slice_length = slice_end - slice_start + 1
        
        comment = '{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(part.oldname, part.oldstart, part.oldend, start, end, slice_start, slice_end)
        start_offset = slice_start - part.oldstart
        end_offset = part.oldend - slice_end
        if part.strand == 1:
            newstart = part.newstart + start_offset
            newend = part.newend - end_offset
        elif part.strand == -1:
            newstart = part.newstart + end_offset
            newend = part.newend - start_offset
        newlength = newend - newstart + 1
        parts.append([part.newname, newstart, newend, newlength, part.strand, part.parttype, comment])
    return parts

def order(start, end):
    if start < end:
        return start, end
    else:
        return end, start

def transfer_genome(new, draft, output, prefix):

    transfer = []
    for scaffold in draft:
        for part in draft[scaffold]:
            if part.parttype == 'haplotype' or part.parttype == 'removed':
                transfer.append(part)
                continue
                
            newstart, newend = order(part.newstart, part.newend)
            parts = get_parts(part.newname, newstart, newend, new)
            
            start = part.oldstart
            for p in parts:
                end = start + p[3] - 1
                transfer.append(Part(part.oldname, start, end, p[0], p[1], p[2], part.strand * p[4], p[5]))
                start = end + 1
    
    for scaffold in new:
        if scaffold.startswith(prefix):
            for part in new[scaffold]:
                if part.parttype == 'active':
                    transfer.append(part)
    
    try:
        with open(output, 'w') as o:
            transfer.sort(key=lambda x: (x.oldname, x.oldstart))
            for part in transfer:
                o.write(repr(part) + "\n")
    except IOError:
        print("Can't open output file", output)
        sys.exit()


def get_args():
    parser = argparse.ArgumentParser(description='''Output new database containing old linkage information on new HaploMerger output
    
        -n new_genome
        -d draft
        -o output
        -p prefix
        ''')

    parser.add_argument('-n', '--new_genome', type=str, required=True)
    parser.add_argument('-d', '--draft', type=str, required=True)
    parser.add_argument('-o', '--output', type=str, required=True)
    parser.add_argument('-p', '--prefix', type=str, required=True)
    return parser.parse_args()

if __name__ == '__main__':
    
    args = get_args()

    new = load_genome(args.new_genome)
    draft = load_genome(args.draft)
    
    transfer_genome(new, draft, args.output, args.prefix)