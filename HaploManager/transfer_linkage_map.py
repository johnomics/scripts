#!/usr/bin/python3 -u

import os
import sys
import pathlib
import argparse
import glob
import sqlite3 as sql
from pprint import pprint
from Bio import SeqIO
from collections import defaultdict
from operator import itemgetter

class Overlap():
    def __init__(self, name, size, start, end, length, strand = 0, my_id = 0, newname=''):
        self.name = name
        self.size = int(size)
        self.start = int(start)+1
        self.end = int(end)
        self.length = int(length)
        self.strand = int(strand)
        self.id = int(my_id)
        self.newname = newname
    
    def __repr__(self):
        return '{}:{}-{} {}bp {} {}'.format(self.name, self.start, self.end, self.length, self.strand, self.newname)

class Part():
    def __init__(self, line, prefix):
        self.scaffold_id, self.scaffold_len, self.sub_scaffold_id, self.new_portion_id, self.active_portion, \
        old_scaffold1_name, old_scaffold1_size, old_scaffold1_id, start1, end1, strand1, len1, \
        old_scaffold2_name, old_scaffold2_size, old_scaffold2_id, start2, end2, strand2, len2, \
        self.connection, self.score, self.sc_Ns, self.tsc_LCs, self.qsc_Ns, self.qsc_LCs, \
        self.active_portion_updated, self.connection_updated, self.active_portion_manual, self.connection_manual \
            = line.rstrip().split('\t')

        if old_scaffold1_name != '0':
            old_scaffold1_name = prefix + old_scaffold1_name
        
        if old_scaffold2_name != '0':
            old_scaffold2_name = prefix + old_scaffold2_name

        self.scaffold1 = Overlap(old_scaffold1_name, old_scaffold1_size, start1, end1, len1, strand1, old_scaffold1_id)
        self.scaffold2 = Overlap(old_scaffold2_name, old_scaffold2_size, start2, end2, len2, strand2, old_scaffold2_id)


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
    
    def collapse_consecutive(self):
        newparts = []
        for i in range(0, len(self.mapparts)):
            mpi = self.mapparts[i]
            if len(newparts) == 0:
                newparts.append(mpi)
            else:
                last = newparts[-1]
                if (last.scaffold==mpi.scaffold and
                    abs(last.end-mpi.start) == 1 and
                    last.chromosome == mpi.chromosome and
                    last.cm == mpi.cm):
                    last.end = mpi.end
                    last.length = last.length + mpi.length
                else:
                    newparts.append(mpi)

        self.mapparts = newparts

    def collapse(self):
        if len(self.mapparts) == 1:
            return
        
        self.collapse_consecutive()
        
        i = 0
        j = 1
        while i < len(self.mapparts)-1:
            mpi = self.mapparts[i]

            if mpi.chromosome == 0 and mpi.cm == -1:
                i+=1
                continue

            parts = []
            j = i+1
            while True:
                if j == len(self.mapparts):
                    break
                mpj = self.mapparts[j]
                
                # Collapse chromosome-only-first blocks (eg mpi is chr1, -1 cm, mpj is chr1, 0.735 cm; assume mpi should be 0.735 also)
                if mpi.chromosome == mpj.chromosome and mpi.cm == -1 and mpj.cm != -1:
                    mpi.cm = mpj.cm

                if (mpj.chromosome == 0 and mpj.cm == -1) or (mpj.chromosome == mpi.chromosome and (mpj.cm == -1 or mpj.cm == mpi.cm)):
                    parts.append(j)
                    j+=1
                    continue
                
                break
            
            while parts and self.mapparts[parts[-1]].chromosome == 0 and self.mapparts[parts[-1]].cm == -1:
                del parts[-1]
            if parts:
                for p in parts:
                    mpi.length = mpi.length + self.mapparts[p].length
                    mpi.end = self.mapparts[p].end
                for p in reversed(parts):
                    del self.mapparts[p]
            i+=1

        return        
        
    def order(self, start, end):
        if start < end:
            return start, end
        else:
            return end, start

    def get_map_slice(self, s_start, s_end, mp):
        s_start, s_end = self.order(s_start, s_end)
        m_start, m_end = self.order(mp.start, mp.end)
        if m_start > s_end or m_end < s_start: # Map block is outside of scaffold block
            return None, None
        slice_start = max(s_start, m_start)
        slice_end = min(s_end, m_end)
        return slice_start, slice_end

    def get_scaffold_range(self, scaffold, mapparts):
        if '_' in scaffold.name:
            gname, gstart, gend, *args = scaffold.name.split('_')
        else:
            gname, gstart, gend = scaffold.name, scaffold.start, scaffold.end
        gstart = int(gstart)
        gend = int(gend)
        if gstart < gend:
            s_start, s_end = scaffold.start, scaffold.end
        else:
            s_start, s_end = gstart-scaffold.start+1, gstart-scaffold.end+1

        partrange = range(0, len(mapparts))
        if scaffold.strand == -1:
            s_start, s_end = s_end, s_start
            partrange = range(len(mapparts)-1, -1, -1)
    
        return s_start, s_end, partrange

    def slice(self, scaffold, new_start):
        slice_parts = []
        s_start, s_end, partrange = self.get_scaffold_range(scaffold, self.mapparts)
        for p in partrange:
            mp = self.mapparts[p]
            slice_start, slice_end = self.get_map_slice(s_start, s_end, mp)
            if slice_start is not None:
                slice_length = slice_end - slice_start + 1
                comment = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(scaffold.name, scaffold.start, scaffold.end, s_start, s_end, mp.start, mp.end, slice_start, slice_end)
                new_mp = MapPart(scaffold.name, new_start, new_start+slice_length-1, slice_length, mp.chromosome, mp.cm, comment)
                slice_parts.append(new_mp)
                new_start += slice_length

        return slice_parts


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


def load_linkage_map(database, prefix, errorarg):
    conn_in, ci = open_input_database(database)
    
    errors = load_errors(errorarg)

    linkage_map = {}
    for chromosome, cm, scaffold, start, end, length, *args in ci.execute('select * from scaffold_map order by scaffold, start'):

        if scaffold in errors and start in errors[scaffold]:
            chromosome = 0
            cm = -1

        scaffold = prefix + scaffold
        if not scaffold in linkage_map:
            linkage_map[scaffold] = MapScaffold()

        linkage_map[scaffold].append(scaffold, start, end, length, chromosome, cm)
    
    for scaffold in linkage_map:
        linkage_map[scaffold].collapse()

    return linkage_map


def open_genome(genome):
    if not os.path.isdir(genome):
        print("Genome argument is not a directory!")
        sys.exit()
    
    result_dir = glob.glob(genome + "/*x.result")
    
    if len(result_dir) != 1:
        print("Got {} result directories, expected 1".format(len(result_dir)))
        sys.exit()

    if not os.path.isdir(result_dir[0]):
        print("Found {} in {} but it is not a directory!".format(result_dir, genome))
        sys.exit()

    return result_dir[0]


def open_output_database(output):
    try:
        conn = sql.connect(output)
        db = conn.cursor()
        db.execute('drop table if exists scaffold_map')
        db.execute('''create table scaffold_map
                     (chromosome integer, cm real, scaffold text, start integer, end integer, length integer, comments text)''')

        cursor = conn.cursor()
        return conn, cursor
    except sql.Error as sqle:
        print(sqle)
        sys.exit(1)


def load_scaffolds(result_dir, prefix):
    
    if os.path.isfile(result_dir + "/hm.new_scaffolds_updated"):
        scaffolds_path = result_dir + "/hm.new_scaffolds_updated"
    else:
        scaffolds_path = result_dir + "/hm.new_scaffolds"

    if not os.path.isfile(scaffolds_path):
        print("Can't find scaffolds file {}".format(scaffolds_path))
        sys.exit()
    
    scaffolds=[]
    try:
        with open(scaffolds_path) as s:
            for line in s:
                if line.startswith('#'):
                    continue

                f = line.rstrip().split('\t')
                if len(f) <= 1:
                    continue
                
                scaffolds.append(Part(line, prefix))
    except IOError:
        print("Failed to load scaffolds file {}".format(scaffolds_path))
        sys.exit()

    return scaffolds

def load_unpaired(result_dir, prefix):
    unpaired_path = result_dir + "/hm.unpaired_updated"
    if not os.path.isfile(unpaired_path):
        print("Can't find unpaired file {}".format(unpaired_path))
        sys.exit()
    
    unpaired = []
    try:
        with open(unpaired_path) as u:
            for line in u:
                if line.startswith('#'):
                    continue
                
                f = line.rstrip().split('\t')
                if len(f) <= 1:
                    continue
                new_name, old_name, size, start, length, fulltype, ns, lcs = f
                size, start, length = int(size), int(start), int(length)
                old_name = prefix + old_name
                if length >= 500:
                    unpaired.append(Overlap(old_name, size, start, start+length, length, newname=new_name))
    except IOError:
        print("Failed to load unpaired file {}".format(unpaired_path))
    return unpaired


def load_genome(genome, prefix):
    result_dir = open_genome(genome)
    
    scaffolds = load_scaffolds(result_dir, prefix)
    
    unpaired = load_unpaired(result_dir, prefix)
    
    return scaffolds, unpaired

def get_linkage(scaffold, linkage_map, new_start):
    if scaffold.name not in linkage_map:
        portion_map = MapScaffold()
        portion_map.append(scaffold.name, scaffold.start, scaffold.end, scaffold.length, 0, -1)
        return portion_map.slice(scaffold, new_start)
    else:
        return linkage_map[scaffold.name].slice(scaffold, new_start)





def make_new_map(linkage_map, scaffolds, new_map, prefix):
    scaffold_num = -1
    cur_part = ''
    scaffold_start = 1
    for part in scaffolds:
        if cur_part != part.scaffold_id:
            cur_part = part.scaffold_id
            scaffold_num += 1
            scaffold_newname = '{}Sc{:07d}'.format(prefix, scaffold_num)
            scaffold_start = 1
        
        scaffold = None
        active_portion = part.active_portion_manual if part.active_portion_manual != "0" else part.active_portion
        scaffold = part.scaffold1 if active_portion == "1" else part.scaffold2
        map_parts = get_linkage(scaffold, linkage_map, scaffold_start)
        
        for mp in map_parts:
            new_map.append([mp.chromosome, mp.cm, scaffold_newname, mp.start, mp.end, mp.length, mp.comment])
        scaffold_start += scaffold.length

def add_unpaired(linkage_map, unpaired, new_map, prefix):
    for part in unpaired:
        map_parts = get_linkage(part, linkage_map, 1)
        for mp in map_parts:
            newname = prefix + part.newname
            new_map.append([mp.chromosome, mp.cm, newname, mp.start, mp.end, mp.length, mp.comment])

def load_additional(linkage_map, additional, prefix):
    if not additional:
        return
    
    try:
        with open(additional, 'r') as f:
            sequences = SeqIO.to_dict(SeqIO.parse(f, "fasta"))
    except IOError:
        print("Can't load additional genome from file {}!".format(fasta))
        sys.exit()

    for scaffold, seq in sequences.items():
        scaffold = prefix + scaffold
        if scaffold not in linkage_map:
            linkage_map[scaffold] = MapScaffold()
            linkage_map[scaffold].append(scaffold, 1, len(seq), len(seq), 0, -1)

def fill_map(mapparts):
    i = 0
    j = 1

    while i < len(mapparts) - 1:
        mpi = mapparts[i]
    
        if mpi[0] == 0 or mpi[1] == -1:
            i += 1
            continue
        
        parts = []
        j = i+1
        while j < len(mapparts):
            mpj = mapparts[j]
            # Different scaffolds
            if mpi[2] != mpj[2]:
                break
            
            # If j chromosome is real and chromosomes or cms don't match, skip
            if mpj[0] != 0 and (mpj[0] != mpi[0] or mpj[1] != mpi[1]):
                break
            
            # If chromosomes and cms match, fill intermediate blocks
            if mpj[0] == mpi[0] and mpj[1] == mpi[1]:
                for p in parts:
                    mapparts[p][0] = mpi[0]
                    mapparts[p][1] = mpi[1]
                break

            parts.append(j)
            j += 1
            
        i += 1


def write_new_map(co, new_map):
    
    new_map.sort(key=itemgetter(2,3))

    fill_map(new_map)

    for mp in new_map:
        co.execute('insert into scaffold_map values (?,?,?,?,?,?,?)', mp)

def get_args():
    parser = argparse.ArgumentParser(description='''Output new database containing old linkage information on new HaploMerger output
    
        -d database
        -g genome
        -a additional
        -o output
        -p prefix
        -e errors
        ''')

    parser.add_argument('-d', '--database', type=str, required=True)
    parser.add_argument('-g', '--genome', type=str, required=True)
    parser.add_argument('-o', '--output', type=str, required=True)
    parser.add_argument('-a', '--additional', type=str, required=False)
    parser.add_argument('-p', '--prefix', type=str, required=True)
    parser.add_argument('-e', '--errors', type=str, required=False)
    return parser.parse_args()

if __name__ == '__main__':
    
    args = get_args()

    linkage_map = load_linkage_map(args.database, args.prefix, args.errors)

    load_additional(linkage_map, args.additional, args.prefix)
    
    scaffolds, unpaired = load_genome(args.genome, args.prefix)
        
    conn_out, co = open_output_database(args.output)

    new_map = []
    make_new_map(linkage_map, scaffolds, new_map, args.prefix)
    
    add_unpaired(linkage_map, unpaired, new_map, args.prefix)

    write_new_map(co, new_map)
    
    conn_out.commit()
    
    conn_out.close()