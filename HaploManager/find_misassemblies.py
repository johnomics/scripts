#!/usr/bin/python3 -u

import os
import sys
import pathlib
import argparse
import re
import string
import sqlite3 as sql
from pprint import pprint
from Bio import SeqIO
from collections import defaultdict, namedtuple
from operator import itemgetter
from termcolor import colored

Gap = namedtuple('Gap', 'scaffold start end')
Marker = namedtuple('Marker', 'maternal paternal intercross')
Gene = namedtuple('Gene', 'start end attributes')
Portion = namedtuple('Portion', 'scaffold start end strand length')

class MapPart:
    def __init__(self, scaffold, start, end, chromosome, cm, parttype, comment=""):
        self.scaffold = scaffold
        self.start = int(start)
        self.end = int(end)
        self.chromosome = chromosome
        self.cm = cm
        self.parttype = parttype
        self.comment = comment

    @property
    def length(self):
        return self.end - self.start + 1
    def __repr__(self):
        return '{}:{:>8}-{:>8} {:>8}bp\t{}\t{}\t{}\t{}'.format(self.scaffold, self.start, self.end, self.length, self.chromosome, self.cm, self.parttype, self.comment)

class Links:
    def __init__(self):
        self.next_cm = -1
        self.prev_cm = -1

class Correction:
    def __init__(self, scaffold, start, end, breaktype, details=""):
        self.scaffold = scaffold
        self.start = int(start)
        self.end = int(end)
        self.breaktype = breaktype
        self.details = details

    def __repr__(self):
        return '{}\t{}\t{}\t{}\t{}'.format(self.scaffold, self.start, self.end, self.breaktype, self.details)
        
def load_corrections(correctfile, scaffolds):
    if not correctfile:
        return []
    corrections = defaultdict(lambda:defaultdict(Correction))
    try:
        with open(correctfile, 'r') as c:
            for line in c:
                mode, scaffold, start, end, breaktype, *args = line.rstrip().split('\t')
                details = ""
                if args:
                    details = args[0]
                if end == 'End':
                    end = scaffolds[scaffold][max(scaffolds[scaffold])].end
                corrections[scaffold][int(start)] = Correction(scaffold, start, end, breaktype, details)
        return corrections
    except IOError:
        print("Can't open corrections file", correctfile)
        sys.exit()

class AGP:
    def __init__(self, scaffold, start, end, part, parttype):
        self.scaffold = scaffold
        self.start = int(start)
        self.end = int(end)
        self.part = part
        self.parttype = parttype
        
    def __repr__(self):
        return '{}\t{}\t{}\t{}\t{}\t{}'.format(self.scaffold, self.start, self.end, self.part, self.parttype, self.length)
        
    @property
    def length(self):
        return self.end-self.start+1

def load_agp(agp):
    
    scaffolds = defaultdict(lambda: defaultdict(AGP))
    broken_scaffolds = {}
    try:
        with open(agp, 'r') as a:
            for line in a:
                scaffold, start, end, part, parttype, *args = line.rstrip().split('\t')
                if ':' in args[-1]:
                    old_scaffold, oldstart, oldend, oldpart = args[-1].split(':')
                    if old_scaffold not in broken_scaffolds:
                        broken_scaffolds[old_scaffold] = {}
                    if scaffold not in broken_scaffolds[old_scaffold]:
                        broken_scaffolds[old_scaffold][scaffold] = 0
                    broken_scaffolds[old_scaffold][scaffold] += 1
                scaffolds[scaffold][int(start)] = AGP(scaffold, start, end, part, parttype)
        
        return scaffolds, broken_scaffolds
    except IOError:
        print("Can't open AGP file", agp)
        sys.exit()

Haplotype = namedtuple("Haplotype", "hapname hapstart hapend hapstrand name start end strand")

class Merge:
    def __init__(self, scaffold, start, end, new, newstart, newend, strand, parttype, hap=None):
        self.scaffold = scaffold
        self.start = int(start)
        self.end = int(end)
        self.new = new
        self.newstart = int(newstart)
        self.newend = int(newend)
        self.strand = strand
        self.parttype = parttype
        self.hap = None
        if hap:
            self.hap = Haplotype(hap[0], hap[1], hap[2], hap[3], None, None, None, None)
        
    def __repr__(self):
        out = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(self.scaffold, self.start, self.end, self.new, self.newstart, self.newend, self.strand, self.parttype)
        if self.hap:
            out += '\t{}\t{}\t{}\t{}'.format(self.hap.hapname, self.hap.hapstart, self.hap.hapend, self.hap.hapstrand)
        return out
        
def load_merged(merged):
    merged_scaffolds = defaultdict(lambda: defaultdict(Merge))
    try:
        with open(merged, 'r') as m:
            for line in m:
                scaffold, start, end, new, newstart, newend, strand, parttype, *args = line.rstrip().split('\t')
                merged_scaffolds[scaffold][int(start)] = Merge(scaffold, start, end, new, newstart, newend, strand, parttype, args)
        return merged_scaffolds
    except IOError:
        print("Can't open merged TSV file", merged)
        sys.exit()

def load_gff(gff):
    genes = defaultdict(list)
    if not gff:
        return genes
    try:
        with open(gff) as g:
            for line in g:
                if line.startswith('#'):
                    continue
                f = line.rstrip().split('\t')
                scaffold, source, featuretype, start, end, score, strand, phase, attributes = f
                if featuretype == 'gene':
                    genes[scaffold].append(Gene(int(start), int(end), attributes))
    except IOError:
        print("Failed to load GFF file {}".format(gff))
    
    return genes


def open_input_database(database):
    try:
        if not os.path.exists(database):
            raise IOError
        conn = sql.connect(database)
        cursor = conn.cursor()
        return conn, cursor
    except IOError:
        print("Can't open database {}".format(database))
        sys.exit(1)
    except sql.Error:
        print("SQL error")
        sys.exit(1)

def load_linkage_map(database):
    conn_in, ci = open_input_database(database)
    
    genome = {}
    linkage_map = {}
    for chromosome, cm, scaffold, start, end, length, parttype, comment in ci.execute('select * from scaffold_map order by scaffold, start'):
        if not scaffold in genome:
            genome[scaffold] = []
        genome[scaffold].append(MapPart(scaffold, start, end, chromosome, cm, parttype, comment))

        if not chromosome in linkage_map:
            linkage_map[chromosome] = {}
        if cm != -1 and not cm in linkage_map[chromosome]:
            linkage_map[chromosome][cm] = Links()

    return genome, linkage_map


def merge_broken(broken_scaffolds, genome, corrections, scaffolds):
    
    broken = 0
    for old in sorted(broken_scaffolds):
        if len(broken_scaffolds[old]) > 1:
            broken += 1
            carryon = False
            for new in broken_scaffolds[old]:
                if new in genome:
                    carryon = True
            if not carryon:
                return
            print(old)
            for new in sorted(broken_scaffolds[old]):
                print("\t", new)
                if new in genome:
                    for p in genome[new]:
                        print("\t\t", p)
                        if new in corrections and p.start in corrections[new]:
                            print(colored("\t\t\t" + repr(corrections[new][p.start]), 'red', attrs=['bold']))
                else:
                    print("\t\tMissing!")
                for start in sorted(scaffolds[new]):
                    print(colored("\t\t\t" + repr(scaffolds[new][start]), 'blue', attrs=['bold']))
            for new in sorted(broken_scaffolds[old]):
                if new in corrections:
                    for start in corrections[new]:
                        print("\t\t", corrections[new][start])


def set_up_links(linkage_map):
    for chromosome in linkage_map:
        cms = sorted(linkage_map[chromosome].keys())
        for i, cm in enumerate(cms):
            if i > 0:
                linkage_map[chromosome][cm].prev_cm = cms[i-1]
            if i < len(cms)-1:
                linkage_map[chromosome][cm].next_cm = cms[i+1]

class Part:
    def __init__(self, scaffold, start, end, active, passive):
        self.scaffold = scaffold
        self.start = start
        self.end = end
        self.active = active
        self.passive = passive

    def __repr__(self):
        active = "{}\t{}\t{}\t{}\t{}".format(self.active.scaffold, self.active.start, self.active.end, self.active.strand, self.active.length)
        passive = "{}\t{}\t{}\t{}\t{}".format(self.passive.scaffold, self.passive.start, self.passive.end, self.passive.strand, self.passive.length)
        return "{}\t{}\t{}\t{}\t{}".format(self.scaffold, self.start, self.end, active, passive)

def load_hm_new_scaffolds(pacbio):
    new_scaffolds = defaultdict(list)
    
    if not pacbio:
        return new_scaffolds
    
    try:
        with open(pacbio, 'r') as pb:
            curid = -1
            scfnum = 0
            part_start = 1
            for line in pb:
                if line.startswith(("#","\n")):
                    continue
                f = line.rstrip().split('\t')
                if curid == -1:
                    curid = f[0]
                if f[0] != curid:
                    scfnum += 1
                    curid = f[0]
                    part_start = 1
                active_portion = f[-2] if f[-2] != '0' else f[4]
                if active_portion not in ['1', '2']:
                    print("Invalid active portion!")
                    print(line)

                portion1 = Portion(f[5], int(f[8]), int(f[9]), int(f[10]), int(f[11]))
                portion2 = Portion(f[12], int(f[15]), int(f[16]), int(f[17]), int(f[18]))
                if active_portion == '1':
                    active, passive = portion1, portion2
                else:
                    active, passive = portion2, portion1

                part_end = part_start + active.length - 1
                new_scaffolds[scfnum].append(Part(scfnum, part_start, part_end, active, passive))
                part_start = part_end + 1

        return new_scaffolds
        
    except IOError:
        print("Cannot open PacBio new scaffolds file!", pacbio)
        sys.exit()

def get_linkage_trio(parts, i):
    trio = []
    part = 1
    while len(trio) < 3 and i+part <= len(parts)-1:
        if parts[i+part].cm != -1 and (not trio or parts[i+part].cm != trio[-1][1]):
            trio.append((parts[i+part].chromosome,parts[i+part].cm))
        part += 1
    if len(trio) < 3:
        trio = []
    return trio
    
def find_linkage_errors(scaffold, genome, linkage_map):
    last_trio = []
    for i in range(0, len(genome[scaffold])-2):
        trio = get_linkage_trio(genome[scaffold], i)
        if not trio:
            return
        if trio[0] == trio[2] and trio[0] != trio[1]:
            if trio == last_trio:
                continue
            print("Linkage error:", scaffold, trio)
            last_trio = trio
    return


def find_scaffold_misassemblies(scaffold, genome, linkage_map, corrections):
    misparts = defaultdict(lambda: defaultdict(MapPart))
    misgroups = []
    i = 0
    while i < len(genome[scaffold])-1:
        pi = genome[scaffold][i]
        if pi.chromosome == 0:
            i += 1
            continue
        
        j = i + 1
        parts = [pi]
        while j < len(genome[scaffold]):
            pj = genome[scaffold][j]
            parts.append(pj)
            if pj.chromosome != 0 and pj.cm != -1:
                misassembly = False
                if pi.chromosome != pj.chromosome:
                    misassembly = True
                elif pi.cm != -1 and pj.cm != -1 and pi.cm != pj.cm:
                    lmi = linkage_map[pi.chromosome][pi.cm]
                    if lmi.prev_cm != pj.cm and lmi.next_cm != pj.cm:
                        misassembly = True
                if misassembly:
                    if (scaffold not in corrections or (scaffold in corrections and
                        parts[1].start in corrections[scaffold] and
                        (corrections[scaffold][parts[1].start].breaktype in ['?', 'L'] or
                         corrections[scaffold][parts[1].start].breaktype == 'R' and ',' in corrections[scaffold][parts[1].start].details))):
                        misassembled = True
                        misgroups.append(parts)
                        for p in parts:
                            misparts[p.scaffold][p.start] = p
                break
            j += 1
        i += 1
    
    return misparts, misgroups

def output_scaffold(parts, misparts):
    for p in parts:
        if p.scaffold in misparts and p.start in misparts[p.scaffold]:
            print(colored(p, 'red', attrs=['bold']))
        else:
            print(p)
    
def output_agp_parts(gap, scaffolds):
    agp_parts = []
    print("M\t{}\t{}\t{}\t".format(gap.scaffold, gap.start, gap.end))
    for start in sorted(scaffolds[gap.scaffold]):
        if start >= gap.start and scaffolds[gap.scaffold][start].end <= gap.end:
            print(scaffolds[gap.scaffold][start])
            agp_parts.append(scaffolds[gap.scaffold][start])
    return agp_parts

def output_merged_parts(gap, merged):
    for start in sorted(merged[gap.scaffold]):
        if gap.start <= start <= gap.end or gap.start <= merged[gap.scaffold][start].end <= gap.end: 
            print(colored(merged[gap.scaffold][start], 'blue', attrs=['bold']))
        else:
            print(merged[gap.scaffold][start])

def output_haplotypes(gap, merged, misparts):
    hapdict = {}
    haps = []
    for start in sorted(misparts[gap.scaffold]):
        scaffold, s1,e1, s2,e2, start, end, strand = misparts[gap.scaffold][start].comment.split('\t')
        if scaffold in merged:
            for start in sorted(merged[scaffold]):
                if merged[scaffold][start].hap:
                    hap = repr(merged[scaffold][start])
                    if hap in hapdict:
                        continue
                    hapdict[hap] = 1
                    haps.append(hap)
    for hap in haps:
        print(hap)


def output_genes(gap, genes):
    gene_found = False
    if gap.scaffold not in genes:
        return

    for g in genes[gap.scaffold]:
        if gap.start <= g.start <= gap.end or gap.start <= g.end <= gap.end:
            print("Gene", g.start, g.end, g.attributes)
            gene_found = True
    if not gene_found:
        print("No genes")

def get_intercross(maternal, paternal):
    intercross = ""
    paternal = ''.join([c for c in paternal if c in ['A', 'B', 'H']]) # Strip colors
    for i in range(len(maternal)):
        if maternal[i] == paternal[i]:
            intercross += paternal[i]
        else:
            intercross += 'H'
    return intercross

def load_chromosome_map(ci):
    markers = defaultdict(lambda: defaultdict(namedtuple))

    for chromosome, maternal, cm, paternal in ci.execute('select distinct chromosome, print, cm, clean from chromosome_map order by chromosome, cm'):
        intercross = get_intercross(maternal, paternal)
        markers[chromosome][cm]=Marker(maternal, paternal, intercross)
        if -1 not in markers[chromosome]:
            empty = " " * len(maternal)
            markers[chromosome][-1]=Marker(maternal, empty, empty)

    return markers
    
def print_marker(part, markers):
    marker = markers[part.chromosome][part.cm]
    print("{}\t{}\t{}\tMaternal".format(part.chromosome, part.cm, marker.maternal))
    print("{}\t{}\t{}\tPaternal".format(part.chromosome, part.cm, marker.paternal))
    print("{}\t{}\t{}\tIntercross".format(part.chromosome, part.cm, marker.intercross))
    print()
    
def output_markers(group, gap, markers, agp_parts, ci):

    if not ci:
        return

    top_marker = markers[group[0].chromosome][group[0].cm]
    bottom_marker = markers[group[-1].chromosome][group[-1].cm]
    print_marker(group[0], markers)
    print_marker(group[-1], markers)
    statement = 'select position, pattern, consensus, marker_type, parent_gt, parent_gqs, parent_dps, mq, fs from markers where scaffold=\"{}\" and position>={} and position<={} order by position'.format(
                    gap.scaffold, gap.start-1, gap.end+1)
    current_position = 0
    for position, pattern, consensus, markertype, parent_gt, parent_gqs, parent_dps, mq, fq in ci.execute(statement):
        
        for part in agp_parts:
            if part.start > current_position and part.end < position:
                print()
                print(part)
                print_marker(group[0], markers)
                print_marker(group[-1], markers)
        
        current_position = position
        
        print('\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(position, pattern, consensus, markertype, parent_gt, parent_gqs, parent_dps, mq, fq))

    print_marker(group[0], markers)
    print_marker(group[-1], markers)

def output_new_scaffolds(g, new_scaffolds):
    
    if 'x' in g.scaffold:
        return

    scaffold_num = -1
    name_match = re.compile(r'(.+)Sc(\d+)').match(g.scaffold)
    if name_match:
        scaffold_num = int(name_match.group(2))

    misparts = []
    if scaffold_num in new_scaffolds:
        for part in new_scaffolds[scaffold_num]:
            if g.start <= part.start <= g.end or g.start <= part.end <= g.end:
                misparts.append(part)

    i = 0
    while i < len(misparts)-1:
        if i <= len(misparts)-2:
            j = i+1
            if misparts[i].active.scaffold != misparts[j].active.scaffold and misparts[i].passive.scaffold == misparts[j].passive.scaffold:
                print(misparts[i])
                print(misparts[j])
                breakpoint = misparts[j].passive.end if misparts[j].passive.strand == -1 else misparts[j].passive.start+1
                print("P\t{}\t{}\t{}".format(misparts[j].passive.scaffold, 'B', breakpoint))
        if i <= len(misparts)-3:
            k = i+2
            if misparts[i].passive.scaffold == misparts[j].active.scaffold == misparts[k].passive.scaffold:
                print(misparts[i])
                print(misparts[j])
                print(misparts[k])
                print("P\t{}\t{}\t{}-{}".format(misparts[j].active.scaffold, 'R', misparts[j].active.start+1, misparts[j].active.end))

        i += 1

def find_misassemblies(genome, linkage_map, scaffolds, merged, genes, corrections, snps, new_scaffolds):
    
    if snps:
        snpdb, snpio = open_input_database(snps)
        markers = load_chromosome_map(snpio)
    else:
        snpio = None
        markers = []

    for scaffold in sorted(genome):
        find_linkage_errors(scaffold, genome, linkage_map)
        misparts, misgroups = find_scaffold_misassemblies(scaffold, genome, linkage_map, corrections)
        
        if misgroups:
            output_scaffold(genome[scaffold], misparts)
            for group in misgroups:
                if group[1].start == group[-2].end + 1:
                    g = Gap(scaffold, group[0].start, group[1].end)
                else:
                    g = Gap(scaffold, group[0].start, group[-1].end)
                agp_parts = output_agp_parts(g, scaffolds)
                output_merged_parts(g, merged)
                output_haplotypes(g, merged, misparts)
                output_genes(g, genes)
                output_markers(group, g, markers, agp_parts, snpio)
                output_new_scaffolds(g, new_scaffolds)
            print("-------")

def get_args():
    parser = argparse.ArgumentParser(description='''Output new database containing old linkage information on new HaploMerger output
    
        -d database
        -s snps
        -a AGP
        -m merged
        -g gff
        -c corrections
        -p pacbio_new_scaffolds
        ''')

    parser.add_argument('-d', '--database', type=str, required=True)
    parser.add_argument('-a', '--agp', type=str, required=True)
    parser.add_argument('-m', '--merged', type=str, required=True)
    parser.add_argument('-c', '--corrections', type=str, required=False)
    parser.add_argument('-g', '--gff', type=str, required=False)
    parser.add_argument('-s', '--snps', type=str, required=False)
    parser.add_argument('-p', '--pacbio', type=str, required=False)

    return parser.parse_args()

if __name__ == '__main__':
    
    args = get_args()

    scaffolds, broken_scaffolds = load_agp(args.agp)

    corrections = load_corrections(args.corrections, scaffolds)

    genes = load_gff(args.gff)
    
    merged = load_merged(args.merged)
    
    genome, linkage_map = load_linkage_map(args.database)

    merge_broken(broken_scaffolds, genome, corrections, scaffolds)
    
    set_up_links(linkage_map)
    
    new_scaffolds = load_hm_new_scaffolds(args.pacbio)
    
    find_misassemblies(genome, linkage_map, scaffolds, merged, genes, corrections, args.snps, new_scaffolds)