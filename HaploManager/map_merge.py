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
    def __init__(self, line, num=-1, prefix = ""):
        self.scaffold_id, self.scaffold_len, self.sub_scaffold_id, self.new_portion_id, self.active_portion, \
        old_scaffold1_name, old_scaffold1_size, old_scaffold1_id, start1, end1, strand1, len1, \
        old_scaffold2_name, old_scaffold2_size, old_scaffold2_id, start2, end2, strand2, len2, \
        self.connection, self.score, self.sc_Ns, self.tsc_LCs, self.qsc_Ns, self.qsc_LCs, \
        self.active_portion_updated, self.connection_updated, self.active_portion_manual, self.connection_manual \
            = line.rstrip().split('\t')

        if num >= 0:
            self.scaffold_name = prefix + 'Sc{:07d}'.format(num)
        self.scaffold1 = Overlap(old_scaffold1_name, old_scaffold1_size, start1, end1, len1, strand1, old_scaffold1_id)
        self.scaffold2 = Overlap(old_scaffold2_name, old_scaffold2_size, start2, end2, len2, strand2, old_scaffold2_id)
        
    def __repr__(self):
        return '{}\t{}'.format(self.scaffold_name, self.scaffold_len)

class OutPart:
    def __init__(self, end, new_name, new_start, new_end, strand, parttype):
        self.end = end
        self.new_name = new_name
        self.new_start = new_start
        self.new_end = new_end
        self.strand = strand
        self.type = parttype
    
    def __repr__(self):
        return '{}\t{}\t{}\t{}\t{}\t{}'.format(self.end, self.new_name, self.new_start, self.new_end, self.strand, self.type)

def load_scaffolds(prefix, new, old_genome, new_genome):
    new_scaffolds_file = "hm.new_scaffolds_updated"
    if not os.path.isfile(new_scaffolds_file):
        new_scaffolds_file = "hm.new_scaffolds"
        if not os.path.isfile(new_scaffolds_file):
            print("Can't find " + new_scaffolds_file)
            sys.exit()
    
    try:
        with open(new_scaffolds_file) as s:
            num = 0
            curid = ''
            curstart = 1
            for line in s:
                if line.startswith('#'):
                    continue

                f = line.rstrip().split('\t')
                if len(f) <= 1:
                    continue
                scaffold_id = f[0]
                if not curid:
                    curid = scaffold_id
                
                if scaffold_id != curid:
                    curid = scaffold_id
                    curstart = 1
                    num += 1

                part = Part(line, num, prefix)

                active_portion = part.active_portion_manual if part.active_portion_manual != '0' else part.active_portion
                scaffold = part.scaffold1 if active_portion == '1' else part.scaffold2
                curend = curstart + scaffold.length - 1
                new_genome[part.scaffold_name][curstart] = OutPart(curend, scaffold.name, scaffold.start, scaffold.end, scaffold.strand, 'merged')
                
                old_genome[scaffold.name][scaffold.start] = OutPart(scaffold.end, part.scaffold_name, curstart, curend, scaffold.strand, "active")
                hapscaffold = part.scaffold1 if active_portion == '2' else part.scaffold2
                if hapscaffold.name != "0":
                    old_genome[hapscaffold.name][hapscaffold.start] = OutPart(hapscaffold.end, part.scaffold_name, curstart, curend, hapscaffold.strand, "haplotype")

                curstart = curend + 1

    except IOError:
        print("Failed to load scaffolds file {}".format(scaffolds_path))
        sys.exit()

def load_refined():
    
    refined_types = defaultdict(str)
    refined_file = "../_hm.unincorpRefiner.log"
    if not os.path.isfile(refined_file):
        print("Can't find _hm.unincorpRefiner.log, no refined types will be output")
        return refined_types
        
    try:
        with open(refined_file) as r:
            for line in r:
                if not (line.startswith('xp') or line.startswith('xf')):
                    continue
                
                new_name, old_name, size, start, length, fulltype, ns, lowcase, ali_len, ali_cov, retained = line.rstrip().split('\t')
                if retained == "1":
                    refined_types[new_name] = "retained"
                else:
                    refined_types[new_name] = "removed"

            return refined_types

    except IOError:
        print("Failed to load refined file {}".format(refined_file))
        sys.exit()

def load_unpaired(prefix, new, old_genome, new_genome):
    
    refined_types = load_refined()
    
    if not os.path.isfile("hm.unpaired_updated"):
        print("Can't find hm.unpaired_updated")
        sys.exit()
    
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

                refined_type = fulltype
                if new_name in refined_types:
                    refined_type = refined_types[new_name]

                new_name = prefix + new_name
                
                overlap = Overlap(old_name, size, start, start+length, length, new_name=new_name)
                new_start = 1
                new_end = overlap.length
                
                new_genome[new_name][new_start] = OutPart(new_end, overlap.name, overlap.start, overlap.end, overlap.strand, refined_type)
                old_genome[overlap.name][overlap.start] = OutPart(overlap.end, new_name, new_start, new_end, overlap.strand, refined_type)
                
    except IOError:
        print("Failed to load unpaired file {}".format(unpaired_path))


def load_merge(prefix, mergedir, old_genome, new_genome):

    results = mergedir + "/genome.genomex.result"
    if not os.path.isdir(results):
        print("Results argument is not a directory!")
        sys.exit()
    
    try:
        os.chdir(results)
    except OSError:
        print("Could not change directory to results directory!")
        sys.exit()


    try:
        new_file = mergedir + "_new.tsv"
        new = open(new_file, 'w')
    except IOError:
        print("Can't open new genome output file!")
        sys.exit()
        
    load_scaffolds(prefix, new, old_genome, new_genome)
    
    load_unpaired(prefix, new, old_genome, new_genome)
    
    return

def collapse(genome):
    for scaffold in genome:
        starts = sorted(genome[scaffold].keys())
        for start_i in starts:
            if start_i not in genome[scaffold]:
                continue
            for start_j in starts:
                if start_j not in genome[scaffold] or start_j <= start_i:
                    continue
                g_i = genome[scaffold][start_i]
                g_j = genome[scaffold][start_j]
                if g_i.new_name == g_j.new_name and g_i.strand == g_j.strand and g_i.type == g_j.type:
                    if g_i.end + 1 == start_j:
                        if g_i.strand == 1 and g_i.new_end + 1 == g_j.new_start:
                            g_i.end = g_j.end
                            g_i.new_end = g_j.new_end
                            del genome[scaffold][start_j]
                        elif g_i.strand == -1 and g_i.new_start == g_j.new_end + 1:
                            g_i.end = g_j.end
                            g_i.new_start = g_j.new_start
                            del genome[scaffold][start_j]
                    
def output_genome(suffix, genome, mergedir):
    
    collapse(genome)
    
    try:
        file = mergedir + "_" + suffix + ".tsv"
        f = open(file, 'w')
    except IOError:
        print("Can't open " + suffix + " genome output file!")
        sys.exit()

    for scaffold in sorted(genome):
        for start in sorted(genome[scaffold]):
            f.write("{}\t{}\t{}\n".format(scaffold, start, genome[scaffold][start]))

def get_args():
    parser = argparse.ArgumentParser(description='''Output map of new scaffolds to old based on HaploMerger results
        -m HaploMerger directory
        -p output prefix
        ''')

    parser.add_argument('-m', '--mergedir', type=str, required=True)
    parser.add_argument('-p', '--prefix', type=str, required=True)
    return parser.parse_args()

if __name__ == '__main__':
    
    args = get_args()

    old_genome = defaultdict(lambda: defaultdict(OutPart))
    new_genome = defaultdict(lambda: defaultdict(OutPart))
    load_merge(args.prefix, args.mergedir, old_genome, new_genome)
    output_genome("old", old_genome, args.mergedir)
    output_genome("new", new_genome, args.mergedir)
