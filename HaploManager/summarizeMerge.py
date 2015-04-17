#!/usr/bin/python3 -u

import os
import sys
import argparse
from collections import defaultdict

class Overlap():
    def __init__(self, name, size, start, end, length, strand = 0, my_id = 0, new_name=''):
        self.name = name
        self.size = int(size)
        self.start = int(start)+1
        self.end = int(end)
        self.length = int(length)
        self.strand = int(strand)
        self.id = int(my_id)
        self.new_name = new_name
    
    def __repr__(self):
        return '{}:{}-{} {}bp {} {}'.format(self.name, self.start, self.end, self.length, self.strand, self.new_name)

class Part():
    def __init__(self, line, num=-1):
        self.scaffold_id, self.scaffold_len, self.sub_scaffold_id, self.new_portion_id, self.active_portion, \
        old_scaffold1_name, old_scaffold1_size, old_scaffold1_id, start1, end1, strand1, len1, \
        old_scaffold2_name, old_scaffold2_size, old_scaffold2_id, start2, end2, strand2, len2, \
        self.connection, self.score, self.sc_Ns, self.tsc_LCs, self.qsc_Ns, self.qsc_LCs, \
        self.active_portion_updated, self.connection_updated, self.active_portion_manual, self.connection_manual \
            = line.rstrip().split('\t')

        if num >= 0:
            self.scaffold_name = 'Sc{:07d}'.format(num)
        self.scaffold1 = Overlap(old_scaffold1_name, old_scaffold1_size, start1, end1, len1, strand1, old_scaffold1_id)
        self.scaffold2 = Overlap(old_scaffold2_name, old_scaffold2_size, start2, end2, len2, strand2, old_scaffold2_id)
        
    def __repr__(self):
        return '{}\t{}'.format(self.scaffold_name, self.scaffold_len)

def load_scaffolds():
    if not os.path.isfile("hm.new_scaffolds"):
        print("Can't find hm.new_scaffolds")
        sys.exit()
    
    scaffolds=[]
    try:
        with open("hm.new_scaffolds") as s:
            num = 0
            curid = ''
            for line in s:
                if line.startswith('#'):
                    continue

                f = line.rstrip().split('\t')
                if len(f) <= 1:
                    continue
                
                scaffolds.append(Part(line, num))
                if not curid:
                    curid = scaffolds[-1].scaffold_id

                if scaffolds[-1].scaffold_id != curid:
                    curid = scaffolds[-1].scaffold_id
                    num += 1

    except IOError:
        print("Failed to load scaffolds file {}".format(scaffolds_path))
        sys.exit()

    return scaffolds

def load_unpaired():
    if not os.path.isfile("hm.unpaired_updated"):
        print("Can't find hm.unpaired_updated")
        sys.exit()
    
    unpaired = []
    try:
        with open("hm.unpaired_updated") as u:
            for line in u:
                if line.startswith('#'):
                    continue
                
                f = line.rstrip().split('\t')
                if len(f) <= 1:
                    continue
                new_name, old_name, size, start, length, fulltype, ns, lcs = f
                size, start, length = int(size), int(start), int(length)
                if length >= 500:
                    unpaired.append(Overlap(old_name, size, start, start+length, length, new_name=new_name))
    except IOError:
        print("Failed to load unpaired file {}".format(unpaired_path))

    return unpaired


def load_genome(results):

    if not os.path.isdir(results):
        print("Results argument is not a directory!")
        sys.exit()
    
    try:
        os.chdir(results)
    except OSError:
        print("Could not change directory to results directory!")
        sys.exit()

    scaffolds = load_scaffolds()
    
    unpaired = load_unpaired()
    
    return scaffolds, unpaired


class OrigPart:
    def __init__(self, start, end, new_name):
        if start > end:
            start, end = end, start
        self.start = start
        self.end = end
        self.length = self.end - self.start + 1
        self.new_name = new_name
    
    def __repr__(self):
        return "{}\t{}\t{}\t{}".format(self.start, self.end, self.length, self.new_name)

def get_original_scaffolds(scaffolds, unpaired):
    original = defaultdict(list)
    for part in scaffolds:
        for scaffold in part.scaffold1, part.scaffold2:
            if scaffold.name == "0":
                continue
            original[scaffold.name].append(OrigPart(scaffold.start, scaffold.end, part.scaffold_name))

    for scaffold in unpaired:
        original[scaffold.name].append(OrigPart(scaffold.start, scaffold.end, scaffold.new_name))

    for scaffold in original:
        original[scaffold] = sorted(original[scaffold], key = lambda x: x.start)

        i = 0
        j = 1
        while j < len(original[scaffold]):
            a = original[scaffold][i]
            merge = True
            while merge:
                merge = False
                b = original[scaffold][j]
                if a.new_name == b.new_name and a.end == b.start - 1:
                    a.end = b.end
                    a.length += b.length
                    original[scaffold].pop(j)
                    merge = False if j == len(original[scaffold]) else True
            i += 1
            j = i+1

    return original

def load_gff(gff):
    genes = defaultdict(list)
    try:
        with open(gff) as g:
            for line in g:
                if line.startswith('#'):
                    continue
                f = line.rstrip().split('\t')
                scaffold, source, featuretype, start, end, score, strand, phase, attributes = f
                if featuretype == 'gene':
                    genes[scaffold].append((int(start), int(end)))
    except IOError:
        print("Failed to load GFF file {}".format(gff))
    
    return genes

def output_breakage_stats(genes, original):
    total_genes = 0
    complete_genes = 0
    for scaffold in genes:
        total_genes += len(genes[scaffold])
        for gene in genes[scaffold]:
            parts = [] # Storing parts and doing double checks so parts of broken genes are available
            for part in original[scaffold]:
                if gene[0] >= part.start and gene[0] <= part.end or gene[1] >= part.start and gene[1] <= part.end:
                    parts.append(part)
            if len(parts) == 1:
                    if gene[0] >= parts[0].start and gene[1] <= parts[0].end:
                        complete_genes += 1
    
    print("Total genes: {}".format(total_genes))
    print("Complete genes: {} ({:5.2f} %)".format(complete_genes, complete_genes/total_genes * 100))
    broken_genes = total_genes - complete_genes
    print("Broken genes: {} ({:5.2f} %)".format(broken_genes, broken_genes/total_genes * 100))

def get_args():
    parser = argparse.ArgumentParser(description='''Output scaffold structure and gene breakages based on HaploMerger merge
    
        -g GFF file
        -r HaploMerger results directory
        ''')

    parser.add_argument('-g', '--gff', type=str, required=True)
    parser.add_argument('-r', '--results', type=str, required=True)
    return parser.parse_args()

if __name__ == '__main__':
    
    args = get_args()

    scaffolds, unpaired = load_genome(args.results)
    
    original = get_original_scaffolds(scaffolds, unpaired)

    genes = load_gff(args.gff)
    
    output_breakage_stats(genes, original)