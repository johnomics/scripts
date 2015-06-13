#!/usr/bin/python3 -u

import os
import sys
import pathlib
import argparse
import glob
import re
import sqlite3 as sql
from pprint import pprint
from Bio import SeqIO
from collections import defaultdict
from operator import itemgetter

class Part():
    def __init__(self, oldname, oldstart, oldend, newname, newstart, newend, strand, parttype, comment="", chromosome="", cm=""):
        self.oldname = oldname
        self.oldstart = int(oldstart)
        self.oldend = int(oldend)
        self.newname = newname
        self.newstart = int(newstart)
        self.newend = int(newend)
        self.strand = int(strand)
        self.parttype = parttype
        self.comment = comment
        self.chromosome = chromosome
        self.cm = cm

    @property
    def length(self):
        return self.oldend - self.oldstart + 1
    
    def __repr__(self):
        return '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(self.oldname, self.oldstart, self.oldend, self.newname,
                                                               self.newstart, self.newend, self.strand, self.parttype,
                                                               self.chromosome, self.cm)

class Comment:
    def __init__(self, name, oldstart, oldend, start, end, slice_start, slice_end, strand):
        self.name = name
        self.oldstart = oldstart
        self.oldend = oldend
        self.start = start
        self.end = end
        self.slice_start = slice_start
        self.slice_end = slice_end
        self.strand = strand

    def __repr__(self):
        return '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(self.name, self.oldstart, self.oldend, self.start, self.end, self.slice_start, self.slice_end, self.strand)

class MapPart:
    def __init__(self, scaffold, start, end, length, chromosome, cm, comment=""):
        self.scaffold = scaffold
        self.start = start
        self.end = end
        self.length = length
        self.chromosome =chromosome
        self.cm = cm
        self.comment = comment
    
    def __repr__(self):
        return '{}:{}-{} {}bp {} {}'.format(self.scaffold, self.start, self.end, self.length, self.chromosome, self.cm)


class MapScaffold:
    def __init__(self):
        self.mapparts = []
    
    def append(self, scaffold, start, end, length, chromosome, cm):
        self.mapparts.append(MapPart(scaffold, start, end, length, chromosome, cm))


def load_genome(mergedgenome):
    genome = defaultdict(list)
    try:
        with open(mergedgenome, 'r') as g:
            for line in g:
                oldname, oldstart, oldend, newname, newstart, newend, strand, typename, *args = line.rstrip().split('\t')
                if typename in ["full", "part"]:
                    typename = "retained"
                genome[oldname].append(Part(oldname, oldstart, oldend, newname, newstart, newend, strand, typename))
    except IOError:
        print("Failed to load genome file {}".format(mergedgenome))
        sys.exit()
    
    return genome

def load_linkage_map(database, errorarg, genome):
    conn_in, ci = open_input_database(database)
    
    errors = load_errors(errorarg)

    linkage_map = {}
    for chromosome, cm, scaffold, start, end, length, *args in ci.execute('select * from scaffold_map order by scaffold, start'):

        if scaffold in errors and start in errors[scaffold]:
            chromosome = 0
            cm = -1

        if not scaffold in linkage_map:
            linkage_map[scaffold] = MapScaffold()

        linkage_map[scaffold].append(scaffold, start, end, length, chromosome, cm)
    
    for scaffold in genome:
        if scaffold not in linkage_map:
            linkage_map[scaffold] = MapScaffold()
            for part in genome[scaffold]:
                if part.parttype in ['active', 'retained']:
                    linkage_map[scaffold].append(scaffold, part.oldstart, part.oldend, part.length, 0, -1)

    return linkage_map, ci

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


def load_errors(errorarg):
    errors = defaultdict(lambda:defaultdict(int))
    if not errorarg:
        return errors

    try:
        if not os.path.isfile(errorarg):
            raise IOError
        with open(errorarg, 'r') as e:
            for line in e:
                scaffold, start = line.rstrip().split('\t')
                start = int(start)
                errors[scaffold][start] = 0

        return errors

    except IOError:
        print("Can't open error file {}".format(errorarg))
        sys.exit(1)

class GFF:
    def __init__(self, line):
        self.scaffold, self.source, self.featuretype, self.start, self.end, self.score, self.strand, self.phase, self.attributes = line.rstrip().split('\t')
        self.start = int(self.start)
        self.end = int(self.end)
    
    def __repr__(self):
        return '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(self.scaffold, self.source, self.featuretype, self.start, self.end, self.score, self.strand, self.phase, self.attributes)

def load_gff(gff):
    genes = []
    try:
        with open(gff) as g:
            for line in g:
                if line.startswith('#') or 'contig' in line:
                    continue
                genes.append(GFF(line))
    except IOError:
        print("Failed to load GFF file {}".format(gff))
    
    return genes

def open_output_database(output):
    try:
        conn = sql.connect(output)
        db = conn.cursor()
        db.execute('drop table if exists scaffold_map')
        db.execute('''create table scaffold_map
                     (chromosome integer, cm real, scaffold text, start integer, end integer, length integer, parttype text, comments text)''')

        cursor = conn.cursor()
        return conn, cursor
    except sql.Error as sqle:
        print(sqle)
        sys.exit(1)


def get_parts(scaffold, start, end, genome):

    parts = []
    haps = []
    for part in genome[scaffold]:
        if part.oldstart > end or part.oldend < start:
            continue
        slice_start = max(part.oldstart, start)
        slice_end = min(part.oldend, end)
        slice_length = slice_end - slice_start + 1
        
        comment = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(part.oldname, part.oldstart, part.oldend, start, end, slice_start, slice_end, part.strand)
        start_offset = slice_start - part.oldstart
        end_offset = part.oldend - slice_end
        if part.strand == 1 or part.strand == 0:
            newstart = part.newstart + start_offset
            newend = part.newend - end_offset
        elif part.strand == -1:
            newstart = part.newstart + end_offset
            newend = part.newend - start_offset
        newpart = Part(part.newname, newstart, newend, part.oldname, part.oldstart, part.oldend, part.strand, part.parttype,
                       Comment(part.oldname, part.oldstart, part.oldend, start, end, slice_start, slice_end, part.strand))
        if part.parttype in ['haplotype', 'gap']:
            haps.append(newpart)
        else:
            parts.append(newpart)
    return parts, haps

def transfer_chromosome_map(ci, co):
    
    table = ci.execute('select * from sqlite_master where name="chromosome_map" and type="table"')
    if len(table.fetchall()) != 1:
        return
    try:
        co.execute('drop table if exists chromosome_map')
        co.execute('''create table chromosome_map
                     (chromosome integer, print text, cm real, original text, clean text, length integer)''')
        
        chrommap = ci.execute('select * from chromosome_map')
        for (chromosome, chromprint, cm, original, clean, length) in chrommap.fetchall():
            co.execute('insert into chromosome_map values (?,?,?,?,?,?)', [chromosome, chromprint, cm, original, clean, length])
    except sql.Error as sqle:
        print(sqle)
        sys.exit(1)

def write_new_map(linkage_map, genome, output, ci):
    conn_out, co = open_output_database(output + "_map.db")

    transfer_chromosome_map(ci, co)

    new_map = []

    genomelength = 0
    for scaffold in linkage_map:
        for part in linkage_map[scaffold].mapparts:
            new_parts, haps = get_parts(scaffold, part.start, part.end, genome)
            for np in new_parts:
                np.chromosome = part.chromosome
                np.cm = part.cm
                new_map.append(np)

    new_map.sort(key=lambda x: (x.oldname, x.oldstart))

    collapse_map(new_map)

    for mp in new_map:
        co.execute('insert into scaffold_map values (?,?,?,?,?,?,?,?)', [mp.chromosome, mp.cm, mp.oldname, mp.oldstart, mp.oldend, mp.oldend-mp.oldstart+1, mp.parttype, repr(mp.comment)])

    conn_out.commit()
    conn_out.close()

def collapse_map(mapparts):
    i = 0
    while i < len(mapparts) - 1:
        mpi = mapparts[i]
        if mpi.chromosome == 0 or mpi.cm == -1:
            i += 1
            continue


        parts = []
        j = i + 1
        merge = False

        while j < len(mapparts):
            mpj = mapparts[j]

            parts.append(j)

            if mpi.chromosome == mpj.chromosome and mpj.cm == -1:
                j += 1
                continue

            # If j chromosome is real and chromosomes or cms don't match, skip
            if mpj.chromosome != 0 and mpj.cm != -1 and  (mpj.chromosome != mpi.chromosome or mpj.cm != mpi.cm):
                if mpi.chromosome == mpj.chromosome: # Even if cms don't match, if chromosomes match, fill them
                    for p in parts:
                        mapparts[p].chromosome = mpi.chromosome
                break

            # If chromosomes and cms match, fill intermediate blocks
            if mpj.chromosome == mpi.chromosome and mpj.cm == mpi.cm:
                for p in parts:
                    mapparts[p].chromosome = mpi.chromosome
                    mapparts[p].cm = mpi.cm
                
                # Split parts into groups belonging to the same scaffolds
                merge_groups = get_merge_groups(mapparts, i, parts)

                # Merge groups belonging to the same scaffolds
                for m in reversed(merge_groups):
                    if len(m) == 1:
                        continue
                    merge=True
                    mpa = mapparts[m[0]]
                    m.pop(0)
                    for mi in m:
                        mpb = mapparts[mi]
                        blockmerge(mpa, mpb)
                    for mi in reversed(m):
                        del mapparts[mi]
                break
            j += 1
        
        if not merge:
            i += 1

    fill_maternal_only_cms(mapparts)
    recollapse(mapparts) # Collapses new filled maternal blocks and blocks with no linkage information

def recollapse(mapparts):
    i = 0
    
    while i < len(mapparts) - 1:
        mpi = mapparts[i]
        j = i + 1
        if j == len(mapparts):
            break
        while (mpi.chromosome == mapparts[j].chromosome and mpi.cm == mapparts[j].cm and
               mpi.oldname == mapparts[j].oldname and mpi.comment.name == mapparts[j].comment.name and
               mpi.oldend == mapparts[j].oldstart-1):
            blockmerge(mpi, mapparts[j])
            del mapparts[j]
        i += 1
    
def fill_maternal_only_cms(mapparts):
    i = 0
    for i in range(len(mapparts)-2):
        for o in ([i, i+1, i+2], [i+2, i+1, i]):
            mp1, mp2, mp3 = mapparts[o[0]], mapparts[o[1]], mapparts[o[2]]
            if mp1.oldname != mp2.oldname or mp1.oldname != mp3.oldname or mp2.oldname != mp3.oldname:
                continue
            if (mp1.chromosome == 0 and mp1.cm == -1 and
                mp2.chromosome != 0 and mp2.cm == -1 and
                mp3.chromosome != 0 and mp3.cm != -1):
                mp2.cm = mp3.cm


def get_merge_groups(mapparts, i, parts):
    merge_groups = [[i]]
    merge_group_i = 0
    current_oldname = mapparts[i].oldname
    current_commentname = mapparts[i].comment.name
    
    for p in parts:
        mapparts[p].chromosome = mapparts[i].chromosome
        mapparts[p].cm = mapparts[i].cm
        if mapparts[p].oldname != current_oldname or mapparts[p].comment.name != current_commentname:
            merge_group_i += 1
            merge_groups.append([])
            current_oldname = mapparts[p].oldname
            current_commentname = mapparts[p].comment.name
        merge_groups[merge_group_i].append(p)
    return merge_groups

def blockmerge(a,b):
    a.oldend = b.oldend
    if a.comment.strand == 1:
        a.comment.oldend = b.comment.oldend
        a.comment.end = b.comment.end
        a.comment.slice_end = b.comment.slice_end
    elif a.comment.strand == -1:
        a.comment.oldstart = b.comment.oldstart
        a.comment.start = b.comment.start
        a.comment.slice_start = b.comment.slice_start

def get_gene_name(feature):
    genename = None
    if feature.featuretype == 'gene':
        genename = re.search(r"ID=(.+?);", feature.attributes).group(1)
    else:
        genename = re.search(r"Parent=(.+?)[;-]", feature.attributes).group(1)
    return genename
    
def set_feature_status(feature, genename, stats, new_parts, haps):
    if len(new_parts) == 1 and not haps:
        if new_parts[0].parttype == 'removed':
            status = 'removed'
        else:
            status = 'active'
    elif not new_parts and len(haps) == 1:
        status = 'haplotype'
    elif len(new_parts) == 1 and len(haps) == 1 and new_parts[0].oldname == haps[0].oldname:
        status = 'overlapped'
    elif not new_parts and not haps:
        status = 'missing'
    else:
        status = 'broken'
    stats[genename][status] += 1
    return genename

def get_gene_status(stats):
    if len(stats) == 1 and 'missing' in stats:
        status = 'missing'
    elif 'broken' in stats or 'overlapped' in stats or 'missing' in stats:
        status = 'broken'
    elif ('haplotype' in stats or 'removed' in stats) and 'active' in stats:
        status = 'hap_broken'
    elif 'haplotype' in stats or 'removed' in stats:
        status = 'removed'
    else:
        status = 'ok'
    return status


def transfer_gff_feature(feature, new_parts, haps):
    part = new_parts[0] if new_parts else haps[0]
    strand = feature.strand
    if part.strand == -1:
        if feature.strand == '+':
            strand = '-'
        elif feature.strand == '-':
            strand = '+'
    output_feature = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(part.oldname, feature.source, feature.featuretype, part.oldstart, part.oldend, feature.score, strand, feature.phase, feature.attributes)
    return output_feature        

def write_new_gff(genes, genome, output):

    stats = defaultdict(lambda: defaultdict(int))
    genefeatures = defaultdict(list)
    genestatus = {}
    
    gfffile = open(output + ".gff", 'w')
    removefile = open(output + '_removed.gff', 'w')
    brokenfile = open(output + '_broken.gff', 'w')
    
    genename = None
    for feature in genes:
        new_parts, haps = get_parts(feature.scaffold, feature.start, feature.end, genome)
        genename = get_gene_name(feature)
        genefeatures[genename].append((feature, new_parts, haps))
        if feature.featuretype != 'gene' and 'RNA' not in feature.featuretype:
            set_feature_status(feature, genename, stats, new_parts, haps)

    genestats = defaultdict(int)
    genetypes = {}
    for gene in stats:

        statuses = sorted(stats[gene].keys())
        genetype = '_'.join(statuses)
        genetypecounts = '_'.join([str(stats[gene][x]) for x in statuses])
        
        if genetype not in genetypes:
            genetypes[genetype] = 0
        genetypes[genetype] += 1
        status = get_gene_status(stats[gene])
        genestats[status] += 1
        
        for feature in genefeatures[gene]:
            if feature[0].featuretype == 'gene':
                feature[0].attributes += 'Note={}:{}'.format(genetype, genetypecounts)
            if status == 'ok':
                output_feature = transfer_gff_feature(feature[0], feature[1], feature[2])
                gfffile.write(output_feature)
            elif status in ['removed', 'missing']:
                removefile.write(repr(feature[0])) # Haplotypes and removed genes
            elif status in ['broken', 'hap_broken']:
                brokenfile.write(repr(feature[0]))
            else:
                print("No status for {}".format(gene))

    gfffile.close()
    removefile.close()
    brokenfile.close()
    
    total = sum(genestats.values())
    print_stat("OK", genestats['ok'], total)
    print_stat("Removed", genestats['removed'], total)
    print_stat("Brokehap", genestats['hap_broken'], total)
    print_stat("Broken", genestats['broken'], total)
    print_stat("Missing", genestats['missing'], total)
    print("Total   genes:\t{}".format(total))

    for genetype in sorted(genetypes, key=genetypes.get, reverse=True):
        print("{}\t{}".format(genetypes[genetype],genetype))

def print_stat(text, stat, total):
    print ('{:<8} genes:\t{:>5}\t{:5.2f} %'.format(text, stat, stat/total*100))

def get_args():
    parser = argparse.ArgumentParser(description='''Output new database containing old linkage information on new HaploMerger output
    
        -m mergedgenome
        -g gff
        -d database
        -e errors
        -o output
        ''')

    parser.add_argument('-m', '--mergedgenome', type=str, required=True)
    parser.add_argument('-g', '--gff', type=str, required=False)
    parser.add_argument('-d', '--database', type=str, required=False)
    parser.add_argument('-e', '--errors', type=str, required=False)
    parser.add_argument('-o', '--output', type=str, required=True)
    return parser.parse_args()

if __name__ == '__main__':
    
    args = get_args()

    genome = load_genome(args.mergedgenome)
    
    if args.database:
        linkage_map, ci = load_linkage_map(args.database, args.errors, genome)
        write_new_map(linkage_map, genome, args.output, ci)

    if args.gff:
        genes = load_gff(args.gff)
        write_new_gff(genes, genome, args.output)