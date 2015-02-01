#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
from collections import defaultdict
from statistics import mean
from os.path import isfile
import sqlite3 as db
from multiprocessing.dummy import Pool as ThreadPool
from itertools import combinations
from copy import deepcopy

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
from Bio.Alphabet import generic_dna
from BCBio import GFF

class GenomeData:
    def __init__(self, args):
        print("Loading genome...")
        overlap_db = self.open_database(args.overlaps)
        overlaps = overlap_db.cursor()
        
        self.sequences = self.load_genome(args.fasta, overlaps)
        
        print("Loading annotation...")
        self.load_annotation(args.gff)
        
        conn = self.open_database(args.database)
        self.db = conn.cursor()
        
        print("Loading blocks...")
        self.blocks = self.load_blocks(self.db, overlaps)
        
        print("Loading errors...")
        self.errors = self.load_errors(args.errors)


    def open_database(self, dbfile):
        try:
            if isfile(dbfile):
                conn = db.connect(dbfile)
            else:
                raise IOError
        except:
            print("Can't open database {}!".format(dbfile))
            sys.exit()
    
        return conn

    def load_genome(self, fasta, overlaps):
        try:
            with open(fasta, 'r') as f:
                sequences = SeqIO.to_dict(SeqIO.parse(f, "fasta"))
        except IOError:
            print("Can't load genome from file {}!".format(fasta))
            sys.exit()

        print("Original stats:", genome_stats([len(sequences[scaffold]) for scaffold in sequences]))

        delete_scaffolds = []
        for scaffold in sequences:
            statement = "select hittype from genome_overlaps where rscaffold=\"{}\"".format(scaffold)
            skip = 0
            for (hittype, ) in overlaps.execute(statement):
                if hittype == '[CONTAINED]':
                    skip = 1
            if skip:
                delete_scaffolds.append(scaffold)
        
        for scaffold in delete_scaffolds:
            del sequences[scaffold]

        print("Remove contained scaffolds:", genome_stats([len(sequences[scaffold]) for scaffold in sequences]))

        return sequences

    def load_annotation(self, gff):
        try:
            with open(gff, 'r') as g:
                for feature in g:
                    if 'gene' in feature:
                        scaffold, source, gfftype, start, end, score, strand, phase, attributes = feature.rstrip().split('\t')
                        if scaffold not in self.sequences:
                            continue
                        featurestrand = None
                        if strand == '+':
                            featurestrand = 1
                        elif strand == '-':
                            featurestrand = -1
                        location = FeatureLocation(ExactPosition(start), ExactPosition(end), strand=featurestrand)
                        feature = SeqFeature(location, type=gfftype)
                        self.sequences[scaffold].features.append(feature)
        except IOError:
            print("Can't load annotation from file {}!".format(gff))
            sys.exit()

    def load_blocks(self, c, overlaps):
        blocks = {}
        block_num = 0
        genome_length = 0

        # Load map blocks from database
        for scaffold, start, end in c.execute("select scaffold, start, end from mapblocks"):
            if scaffold not in self.sequences:
                continue
            if end < start:
                continue
            block_num += 1
            genome_length += end - start + 1
            if scaffold not in blocks:
                blocks[scaffold] = {}
            blocks[scaffold][start] = Block(scaffold, start, end)

        # Store previous and next blocks
        for scaffold in blocks:
            positions = sorted(blocks[scaffold].keys())
            for i in range(0, len(positions)-1):
                if i > 0:
                    blocks[scaffold][positions[i]].prev_block = blocks[scaffold][positions[i-1]].start
                if i < len(positions)-1:
                    blocks[scaffold][positions[i]].next_block = blocks[scaffold][positions[i+1]].start

        print("Loaded {} blocks, length {}".format(block_num, genome_length))
        return blocks

    def load_errors(self, errorfile):
        errors = {}
        if not errorfile:
            return errors
    
        try:
            if isfile(errorfile):
                with open(errorfile) as err:
                    for line in err:
                        scaffold, start = line.rstrip().split('\t')
                        errors[(scaffold, int(start))] = -1
            else:
                raise IOError
        except:
            print("Can't open errors file!")
            sys.exit()
    
        return errors

class PoolStat:
    def __init__(self, count=0, length=0):
        self.count = count
        self.length = length
    
    def __add__(self, other):
        self.count += other.count
        self.length += other.length
        return self

class Stats:
    def __init__(self, name):
        self.name = name
        self.pool_num = 0
        self.log_num = 0
        self.pool_stats = {'orient':PoolStat(), 'order':PoolStat(), 'overlap':PoolStat(), 'error':PoolStat(), 'ok':PoolStat()}
        self.log_count = {'orient':0, 'order':0, 'overlap':0, 'error':0, 'ok':0}
        self.raft_count = {'orient':0, 'order':0, 'overlap':0, 'error':0, 'ok':0}

    def __add__(self, other):
        self.pool_num += other.pool_num
        self.log_num += other.log_num
        for gt in self.pool_stats:
            self.pool_stats[gt] += other.pool_stats[gt]
            self.log_count[gt] += other.raft_count[gt]
            self.raft_count[gt] += other.raft_count[gt]
        return self

    def __repr__(self):
        output = "{}: {:6d} blocks in {:4d} groups\n".format(self.name, self.log_num, self.pool_num)
        output += "Type\tGroups\tLists\tBlocks\tLength\n"
        for gt in 'ok', 'orient', 'order', 'overlap', 'error':
            output += "{}\t{:6d}\t{:6d}\t{:6d}\t{:10d}\n".format(
                  gt, self.pool_stats[gt].count, self.raft_count[gt],
                  self.log_count[gt], self.pool_stats[gt].length)
        return output


class Chromosome:
    def __init__(self, name, genome):
        self.name = str(name)
        self.markers = {}
        self.set_markers(genome.db)
        self.pools = []
        self.genome = genome
        self.mapped_blocks, self.mapped_blocks_length, self.placed_blocks, self.placed_blocks_length = self.set_blocks()
        
    def __repr__(self):
        return('\n'.join([self.name, repr(self.pools), repr(self.markers)]))

    def __iter__(self):
        return iter(self.pools)

    @property
    def stats(self):
        stats = Stats(self.name)
        stats.pool_num = len(self.pools)
        for pool in self.pools:
            gt = pool.pooltype
            for raft in pool:
                stats.log_num += len(raft.logs)
                stats.pool_stats[gt].length += raft.length
                stats.log_count[gt] += len(raft.logs)
                stats.raft_count[gt] += 1
            stats.pool_stats[gt].count += 1
        
        return stats

    def set_markers(self, db):
        prev_cm = -1
        next_cm = -1
        for cm, in db.execute("select distinct cm from scaffold_map where chromosome={} order by cm".format(self.name)):

            if cm == -1:
                self.add_marker(-1, -1, -1)
                continue

            self.add_marker(cm, prev_cm=prev_cm)
            if prev_cm != -1:
                self.update_marker(prev_cm, next_cm = cm)
        
            prev_cm = cm

    def set_blocks(self):
        mapped_blocks = mapped_blocks_length = placed_blocks = placed_blocks_length = 0

        for cm in sorted(self.markers) + [-1]:
            cm_blocks = Pool(self, self.genome)
            cm_block_id = 1
            statement = "select scaffold, start, end, length from scaffold_map where chromosome={} and cm={} order by scaffold, start, end".format(self.name, cm)
            for scaffold, start, end, length in self.genome.db.execute(statement):
                if scaffold not in self.genome.sequences:
                    continue
                if (scaffold, start) in self.genome.errors:
                    continue
                mapped_blocks += 1
                mapped_blocks_length += length
                self.genome.blocks[scaffold][start].add_marker(self.name,cm)

                if cm != -1:
                    placed_blocks += 1
                    placed_blocks_length += length
                cm_blocks.add(Raft(cm_block_id, scaffold, start, self))
                cm_block_id += 1
            if cm != -1:
                self.pools.append(cm_blocks)

        return mapped_blocks, mapped_blocks_length, placed_blocks, placed_blocks_length

    @property
    def marker_chains(self):
        
        marker_chains = {}
        
        for pool in chromosome:
            if pool.marker_chain not in marker_chains:
                marker_chains[pool.marker_chain] = {}
                marker_chains[pool.marker_chain]['count'] = 0
                marker_chains[pool.marker_chain]['type'] = ''

            marker_chains[pool.marker_chain]['type'] = ''
            marker_chains[pool.marker_chain]['count'] += 1
        
        return marker_chains

    def add_marker(self, cm, prev_cm=-1, next_cm=-1):
        self.markers[cm] = Marker(cm, prev_cm, next_cm)
        
    def update_marker(self, cm, prev_cm = -1, next_cm = -1):
        if cm not in self.markers:
            self.add_marker(cm, prev_cm, next_cm)

        self.markers[cm].update_previous(prev_cm)
        self.markers[cm].update_next(next_cm)

    def get_scaffolds(self):
        scaffolds = []
        for pool in self:
            for raft in pool:
                scaffoldlist = []
                
                for scaffold, start, direction in raft.logs:
                    scaffoldlist.append(self.genome.blocks[scaffold][start])
                scaffolds.append(scaffoldlist)
        
        return scaffolds

    def assemble(self, genome, threads):
        for pool in self.pools:
            pool.assemble(pool, cm_merge)
        self.connect(cm_merge)
        for pool in self.pools:
            pool.extend()

        print(self.stats)

#        extend_ok_groups(self, overlaps)


    def connect(self, merger):
        p = 0
        while p < len(self.pools)-1:
            self.pools[p].assemble(self.pools[p+1], merger)
            p = self.split(p)
            if self.pools[p+1]:
                p += 1
            else:
                del self.pools[p+1]


    def split(self, p):
        ordered_rafts = [raft for raft in self.pools[p] if raft.ordered]
        if not ordered_rafts or len(ordered_rafts) == len(self.pools[p]):
            return p

        self.pools.insert(p+1, Pool(self, genome))
        for raft in ordered_rafts:
            self.pools[p+1].add(raft)
            self.pools[p].remove(raft)
        return p+1


class Pool:
    def __init__(self, chromosome, genome = ''):
        self.chromosome = chromosome
        self.genome = genome
        self.rafts = set()

    def __repr__(self):
        output = 'Type: {}\n'.format(self.pooltype)
        for raft in self:
            output += repr(raft) + '\n'
            output += 'Length: {}\n'.format(raft.length)
            output += '-----\n'
        output += '====='
        
        return output
    
    def __iter__(self):
        return iter(self.rafts)
    
    def __len__(self):
        return len(self.rafts)

    def add(self, raft):
        self.rafts.add(raft)

    def remove(self, raft):
        self.rafts.remove(raft)
    
    @property
    def marker_chain(self):
        for raft in self.rafts:
            return raft.marker_chain
        return ''
    
    @property
    def pooltype(self):
        if len(self.marker_chain) == 1:
            if len(self.rafts) > 1:
                return 'order'
            else:
                return 'orient'
        else:
            if len(self.rafts) > 1:
                return 'overlap'
            else:
                return 'ok'

    def cleanup(self):
        for raft in [raft for raft in self.rafts if not raft.logs]:
            self.rafts.remove(raft)

    def assemble(self, other, merger):
        for a in self.rafts:
            if not a:
                continue

            repeat = True
            while repeat:
                repeat = False
                for b in other.rafts:
                    if not b or repr(a) == repr(b):
                        continue

                    merge = merger(a, b)
                    if merge:
                        repeat = True
                        a.replace(merge)
                        b.empty()
                        break
        other.cleanup()

    def extend(self):
        for raft in self.rafts:
            raft.extend()

def orient(a, b):
    edge_blocks = [a.start, a.last, b.start, b.last]
    if a.start < b.start:
        a_far, a_near, b_near, b_far = sorted(edge_blocks)
    else:
        b_far, b_near, a_near, a_far = sorted(edge_blocks)
    
    a_forward = a.start <= a.last
    b_forward = b.start <= b.last
    return a_near, b_near, a_forward, b_forward

def merge_logs(this, other):
    if not this:
        return deepcopy(other)

    new = deepcopy(this) + deepcopy(other)
    new.sort(key=lambda x: x[1]) # Sort by start position

    this_dir = this[0][2]
    new = [(scaffold, start, this_dir) for scaffold, start, direction in new]
    if this_dir == -1:
        new.reverse()

    return new

def cm_merge(a, b):
    if a.scaffold != b.scaffold:
        return None
    bridge = []
    scaffold = a.genome.blocks[a.scaffold]

    start, target, a_forward, b_forward = orient(a,b)
    this_cm = scaffold[start].cm
    direction = 1 if start < target else -1

    while True:
        if direction == 1:
            start = scaffold[start].next_block
        else:
            start = scaffold[start].prev_block
        
        if start == target:
            bridge = merge_logs(bridge, a.logs)
            bridge = merge_logs(bridge, b.logs)
            return bridge
            
        if start == 0:
            break
        
        block_i_chr = scaffold[start].chromosome
        block_i_cm = scaffold[start].cm
        if block_i_chr == '0' or block_i_cm == -1 and block_i_chr == a.chromosome.name:
            bridge.append((a.scaffold, start, direction))
            continue
    
        if block_i_chr != a.chromosome.name: # This check comes after the cM check because chr==0 is OK
            break
    
        if block_i_cm == this_cm:
            bridge.append((a.scaffold, start, direction))
    
        if block_i_cm == a.chromosome.markers[this_cm].next_cm:
            this_cm = a.chromosome.markers[this_cm].next_cm
            bridge.append((a.scaffold, start, direction))
    
        
        if block_i_cm == a.chromosome.markers[this_cm].prev_cm:
            this_cm = a.chromosome.markers[this_cm].prev_cm
            bridge.append((a.scaffold, start, direction))
        
        if bridge and bridge[-1][1] == start:
            continue

        break

    return None

class Raft:
    def __init__(self, bl_id, scaffold, start, chromosome):
        self.chromosome = chromosome
        self.genome = chromosome.genome
        self.id = bl_id
        self.logs = []
        self.manifest = []
        self.append(scaffold, start, 1)

    def __repr__(self):
        return '\n'.join([repr(m) for m in self.manifest])

    def __hash__(self):
        return self.id

    def __eq__(self, other):
        return self.id == other.id

    def __iter__(self):
        return iter(self.logs)
    
    def __len__(self):
        return len(self.logs)
    
    def __contains__(self, start):
        for log in self.logs:
            if log[1] == start:
                return True
        return False

    @property
    def scaffold(self):
        scaffolds = {m.scaffold:0 for m in self.manifest}
        return '_'.join(sorted(scaffolds.keys()))
        
    @property
    def name(self):
        summarydict = defaultdict(list)
        scaffolds = []
        for m in self.manifest:
            if m.scaffold not in summarydict:
                scaffolds.append(m.scaffold)
            summarydict[m.scaffold].append((m.start,m.end))
        names = []
        for scaffold in scaffolds:
            names.append('{}_{}_{}'.format(scaffold, summarydict[scaffold][0][0], summarydict[scaffold][-1][-1]))
        name = '-'.join(names)
        return name

    @property
    def start(self):
        return self.manifest[0].start
    
    @property
    def last(self):
        return self.logs[-1][1]

    @property
    def end(self):
        return self.manifest[-1].end

    @property
    def length(self):
        length = 0
        for m in self.manifest:
            length += m.length
        return length
    
    @property
    def sequence(self):
        return sum([m.sequence for m in self.manifest], Seq("", generic_dna))

    @property
    def marker_chain(self):
        marker_chain = []
        for m in self.manifest:
            if m.cm != -1:
                marker_chain.append(m.cm)
                if len(marker_chain)>1 and marker_chain[-2] == marker_chain[-1]:
                    del marker_chain[-1]

        if not self.check_chain(marker_chain):
            print("To fix:", marker_chain, "\n", self, "\n")

        return tuple(marker_chain)

    @property
    def ordered(self):
        return len(self.marker_chain) > 1


    def empty(self):
        self.logs = []
        self.manifest = []
        self.update()
        
    def replace(self, logs):
        self.empty()
        for scaffold, start, direction in logs:
            self.append(scaffold, start, direction)
        self.update()

    def append(self, scaffold, start, direction=1):
        self.logs = self.logs + [(scaffold, start, direction)]
        self.manifest = self.manifest + [SummaryBlock(scaffold, start, self.genome.blocks[scaffold][start], direction)]
        self.update()

    def prepend(self, scaffold, start, direction=1):
        self.logs = [(scaffold, start, direction)] + self.logs
        self.manifest = [SummaryBlock(scaffold, start, self.genome.blocks[scaffold][start], direction)] + self.manifest
        self.update()

    def merge(self, other):
        for scaffold, start, direction in other.logs:
            self.append(scaffold, start, direction)

    def update(self):
        self.collapse()

    def check_chain(self, chain):
        
        if len(chain) == 1:
            return True
        
        if chain[0] < chain[1]:
            direction = 1
        else:
            direction = -1
        
        for i in range(1, len(chain)-1):
            if chain[i] < chain[i+1]:
                this_dir = 1
            else:
                this_dir = -1
            if direction != this_dir:
                return False
                break
        
        return True


    def collapse_consecutive(self):
        newsummary = []
        for i in range(0, len(self.manifest)):
            sbi = self.manifest[i]
            if len(newsummary) == 0:
                newsummary.append(sbi)
            else:
                last = newsummary[-1]
                if (last.scaffold==sbi.scaffold and
                    abs(last.end-sbi.start) == 1 and
                    last.cm == sbi.cm):
                    last.end = sbi.end
                else:
                    newsummary.append(sbi)

        self.manifest = newsummary

    def collapse_trios(self):
        if len(self.manifest) < 3:
            return

        newsummary = []
        i = 0
        while i < len(self.manifest):
            sbi = self.manifest[i]
            if i < len(self.manifest)-2:
                sbj = self.manifest[i+1]
                sbk = self.manifest[i+2]
                if (sbi.scaffold == sbj.scaffold and sbi.scaffold == sbk.scaffold and
                    abs(sbi.end - sbj.start) == 1 and abs(sbj.end - sbk.start) == 1 and
                    sbi.cm == sbk.cm             and sbj.cm == -1):
                        sbi.end = sbk.end
                        newsummary.append(sbi)
                        i += 3
                        continue

            newsummary.append(sbi)
            i += 1
        
        self.manifest = newsummary
        

    def collapse(self):
        self.collapse_consecutive()
        self.collapse_trios()

    def order(self):
        if self.manifest[0].cm != -1 and self.manifest[-1].cm != -1 and self.manifest[0].cm > self.manifest[-1].cm:
            self.logs.reverse()
            for i in range(0, len(self.logs)):
                (scaffold, start, direction) = self.logs[i]
                self.logs[i] = (scaffold, start, -1)

            self.manifest.reverse()
            for sb in self.manifest:
                sb.start, sb.end = sb.end, sb.start


    def extend(self):
        self.extend_dir(0)
        self.extend_dir(-1)

    def extend_dir(self, item):
        first_scaffold, first_start, direction = self.logs[item]
        ext_start = first_start
        extend = 0
        starts_to_extend = []
        scaffold = self.genome.blocks[first_scaffold]
        chromosome = scaffold[first_start].chromosome
        
        while True:
            if item == 0 and direction == 1 or item == -1 and direction == -1:
                ext_start = scaffold[ext_start].prev_block
            else:
                ext_start = scaffold[ext_start].next_block
        
            # Extending to ends is OK
            if ext_start == 0:
                break
        
            # If we reach another cM, abandon extension
            if scaffold[ext_start].cm != -1 or scaffold[ext_start].chromosome != '0' and scaffold[ext_start].chromosome != chromosome:
                starts_to_extend = []
                break
        
            starts_to_extend.append(ext_start)
        
        if starts_to_extend:
            if item == 0:
                for start in starts_to_extend:
                    self.prepend(first_scaffold, start, direction)
            else:
                for start in starts_to_extend:
                    self.append(first_scaffold, start, direction)

class SummaryBlock:
    def __init__(self, scaffold, start, block, direction):
        self.scaffold = scaffold
        if direction == 1:
            self.start = start
            self.end = block.end
        elif direction == -1:
            self.start = block.end
            self.end = start
        self.cm = block.cm

    @property
    def length(self):
        return max(self.start,self.end)-min(self.start,self.end)+1

    @property
    def sequence(self):
        if self.start < self.end:
            seq = genome.sequences[self.scaffold][self.start-1:self.end]
            seq.id = '_'.join([self.scaffold, str(self.start), str(self.end)])
        else:
            seq = genome.sequences[self.scaffold][self.end-1:self.start][::-1]
            seq.id = '_'.join([self.scaffold, str(self.start), str(self.end)])
        seq.name = seq.id
        seq.description = seq.id
        return seq

        
    def __repr__(self):
        return '{}:{}-{} ({}, {} bp)'.format(self.scaffold, str(self.start), str(self.end), str(self.cm), str(self.length))


class Block:
    def __init__(self, scaffold, start, end, prev_block=0, next_block=0, chromosome=0, cm=-1):
        self.scaffold = scaffold
        self.start = start
        self.end = end
        self.prev_block = prev_block
        self.next_block = next_block
        self.chromosome = str(chromosome)
        self.cm = cm
    
    @property
    def length(self):
        return self.end - self.start + 1
    
    def __repr__(self):
        return "{}:{}-{} ({}, {}) [{}, {}]".format(self.scaffold, self.start, self.end, self.prev_block, self.next_block, self.chromosome, self.cm)

    def add_marker(self, chromosome, cm):
        self.chromosome = str(chromosome)
        self.cm = cm


class Marker:
    def __init__(self, cm, prev_cm=-1, next_cm=-1):
        self.cm = cm
        self.prev_cm = prev_cm
        self.next_cm = next_cm

    def __repr__(self):
        return('{}-({},{})'.format(self.cm, self.prev_cm, self.next_cm))

    def update_previous(self, prev_cm):
        if prev_cm != -1:
            self.prev_cm = prev_cm
    
    def update_next(self, next_cm):
        if next_cm != -1:
            self.next_cm = next_cm





def get_args():
    parser = argparse.ArgumentParser(description='''Output new FASTA file based on linkage map.
    
        -d database
        -f FASTA file
        -g GFF file
        -e errors file
        -t threads
        -o overlaps''')

    parser.add_argument('-d', '--database', type=str, required=True)
    parser.add_argument('-f', '--fasta', type=str, required=True)
    parser.add_argument('-g', '--gff', type=str, required=True)
    parser.add_argument('-o', '--overlaps', type=str, required=True)
    parser.add_argument('-e', '--errors', type=str)
    parser.add_argument('-t', '--threads', type=int, default=1)

    return parser.parse_args()

def genome_stats(lengths):
    scaffolds=len(lengths)
    genome_length = sum(lengths)
    lengths.sort(reverse=True)
    N50_threshold = genome_length / 2
    
    N50_sum = 0
    N50 = 0
    for length in lengths:
        N50_sum += length
        if N50_threshold < N50_sum:
            N50 = length
            break
    
    return "Scaffolds: {}\tLength: {}\tN50: {}".format(scaffolds, genome_length, N50)


def add_overlap(rscaffold, rstart, rend, qscaffold, pcid, hittype):
    if rstart < rend:
        location = FeatureLocation(ExactPosition(rstart), ExactPosition(rend), strand=1)
    else:
        location = FeatureLocation(ExactPosition(rend), ExactPosition(rstart), strand=1)
    feature = SeqFeature(location, type="Overlap")
    genome.sequences[rscaffold].features.append(feature)


# Unused; will load from show-coords dump of nucmer delta file
def load_overlaps(overlaps):
    try:
        with open(overlaps, 'r') as o:
            for line in o:
                f = line.rstrip().split('\t')
                if len(f) < 11:
                    continue
                rstart, rend, qstart, qend, rhitlen, qhitlen, pcid, rseqlen, qseqlen, rscaffold, qscaffold, *rem = f
                if rscaffold == qscaffold:
                    continue
                hittype = ''
                if len(rem) == 1:
                    hittype = rem[0]
                
                add_overlap(rscaffold, rstart, rend, qscaffold, pcid, hittype)
                add_overlap(qscaffold, qstart, qend, rscaffold, pcid, hittype)
    except IOError:
        print("Can't load overlaps from file {}!".format(overlaps))
        sys.exit()


def get_genome_overlaps(scaffolds, overlaps):
    found = False
    statement = '''select * from genome_overlaps
                   where rscaffold in ({0}) and qscaffold in ({0})
                     and rscaffold != qscaffold'''.format(
                ','.join('?'*len(scaffolds)))
    for rstart, rend, qstart, qend, rhitlen, qhitlen, pcid, rseqlen, qseqlen, rpccov, qpccov, \
        rscaffold, qscaffold, hittype in overlaps.execute(statement, list(scaffolds)*2):
        found = True
        print(rscaffold, len(genome.sequences[rscaffold]), rstart, rend, rhitlen, rpccov,
              qscaffold, len(genome.sequences[qscaffold]), qstart, qend, qhitlen, qpccov,
              pcid, hittype, sep='\t')
    
    if not found:
        print("No genome overlaps found")


def order_hit_ends(start, end):
    if start > end:
        start, end = end, start
    return start, end


def get_pacbio_hits(scaffolds, overlaps):
    pacbio_scaffolds = {}
    pacbio_lengths = {}
    statement = 'select * from pacbio_overlaps where qscaffold in ({})'.format(','.join('?'*len(scaffolds)))
    for rstart, rend, qstart, qend, rhitlen, qhitlen, pcid, rseqlen, qseqlen, rpccov, qpccov, rscaffold, qscaffold, hittype \
      in overlaps.execute(statement, list(scaffolds)):
        gstart, gend = order_hit_ends(qstart, qend)
        hits_to_keep = {}
        for start in scaffolds[qscaffold]:
            bstart, bend = order_hit_ends(start, scaffolds[qscaffold][start])
#            if (gstart > bstart and gend < bend
#               or gstart > bend or gend < bstart):
            if (gstart > bend or gend < bstart):
                continue
            hits_to_keep[gstart] = gend

        if hits_to_keep:
            if rscaffold not in pacbio_scaffolds:
                pacbio_scaffolds[rscaffold] = defaultdict(list)
                pacbio_lengths[rscaffold] = rseqlen
            pacbio_scaffolds[rscaffold][qscaffold].append((qstart, qend, qhitlen, rstart, rend, pcid))

    return pacbio_scaffolds, pacbio_lengths

def order_pacbio_hits(pbscaffold_hits):
    pacbio_order = {}
    for gs in pbscaffold_hits:
        for (qstart, qend, qhitlen, rstart, rend, pcid) in sorted(pbscaffold_hits[gs], key=lambda x:x[2], reverse=True):
            if rstart > rend:
                qstart, qend = qend, qstart
                rstart, rend = rend, rstart
            ignore = False
            for start, hit in pacbio_order.items():
                if gs == hit['scaffold'] and qstart >= hit['gstart'] and qend <= hit['gend']:
                    ignore = True
            if not ignore:
                pacbio_order[rstart] = {'scaffold': gs, 'pbstart': rstart, 'pbend': rend, 'gstart': qstart, 'gend': qend, 'pcid':pcid}
    return pacbio_order

class Bridge:
    
    def __init__(self, pbs, pblen, hit1, hit2):
        self.pbs = pbs
        self.pblen = pblen
        self.hit1 = hit1
        self.hit2 = hit2
        self.gap = self.hit2['pbstart']-1 - self.hit1['pbend']+1

    def __repr__(self):
        out = ('''{0} ({1:6d} bp) {2} - {3} joins
       {4:16s} ({5:7d} bp) {6:7d} - {7:7d} ({2:7d} - {8:7d}; {9:5.2f}%) and
       {10:16s} ({11:7d} bp) {12:7d} - {13:7d} ({14:7d} - {3:7d}; {15:5.2f}%)'''.format(
                   self.pbs, self.pblen, self.hit1['pbstart'], self.hit2['pbend'],
                   self.hit1['scaffold'], len(genome.sequences[self.hit1['scaffold']]),
                   self.hit1['gstart'], self.hit1['gend'], self.hit1['pbend'], self.hit1['pcid'],
                   self.hit2['scaffold'], len(genome.sequences[self.hit2['scaffold']]),
                   self.hit2['gstart'], self.hit2['gend'], self.hit2['pbstart'], self.hit2['pcid']
               ))
            
        if self.hit1['pbend'] < self.hit2['pbstart']:
            out += "\tFill gap with {} {}-{}".format(self.pbs, self.hit1['pbend']+1, self.hit2['pbstart']-1)
        elif self.hit1['pbstart'] < self.hit2['pbstart'] and self.hit1['pbend'] > self.hit2['pbend']:
            out += "\tSecond hit contained in first hit"
        elif self.hit2['pbstart'] < self.hit1['pbstart'] and self.hit2['pbend'] > self.hit1['pbend']:
            out += "\tFirst hit contained in second hit"
        else:
            out += "\tOverlap of {} bases".format(self.hit1['pbend'] - self.hit2['pbstart'] + 1)
        return out

def add_bridge(bridges, pbs, pblen, hit1, hit2):
    for hit in hit1, hit2:
        if hit['scaffold'] not in bridges:
            bridges[hit['scaffold']] = defaultdict(list)
    bridges[hit1['scaffold']][hit2['scaffold']].append(Bridge(pbs, pblen, hit1, hit2))
    bridges[hit2['scaffold']][hit1['scaffold']].append(Bridge(pbs, pblen, hit2, hit1))
    
def get_pacbio_bridges(scaffolds, pacbio_scaffolds, pacbio_lengths):
    bridges = {}
    for pbs in pacbio_scaffolds:
        if len(pacbio_scaffolds[pbs]) <= 1:
            continue
        pacbio_order = order_pacbio_hits(pacbio_scaffolds[pbs])

        for s1, hit1 in sorted(pacbio_order.items()):
            for s2, hit2 in sorted(pacbio_order.items()):
                if s2 <= s1:
                    continue
                add_bridge(bridges, pbs, pacbio_lengths[pbs], hit1, hit2)

    return bridges


def get_pacbio_overlaps(scaffolds, overlaps):
    
    bridges = {}

    pacbio_scaffolds, pacbio_lengths = get_pacbio_hits(scaffolds, overlaps)
    
    if not pacbio_scaffolds:
        print("No PacBio overlaps found")
        return bridges

    bridges = get_pacbio_bridges(scaffolds, pacbio_scaffolds, pacbio_lengths)

    return bridges

def get_pool_scaffolds(scaffolds, pool, cms=None, fillcms=False):
    for raft in pool:
        for sb in raft.manifest:
            if fillcms and sb.cm != -1 and sb.cm not in cms:
                cms[sb.cm] = 0
            if sb.cm in cms:
                if sb.scaffold not in scaffolds:
                    scaffolds[sb.scaffold] = {}
                if raft.start not in scaffolds[sb.scaffold]:
                    scaffolds[sb.scaffold][raft.start] = raft.end

def get_neighbour(group, bridges, direction):
    neighbour = None
    (thisgroup, ) = group
    if direction == 1:
        this = thisgroup.manifest[-1]
    else:
        this = thisgroup.manifest[0]
    if this.scaffold in bridges:
        for a in bridges[this.scaffold]:
            if this.scaffold != a:
                print(this.scaffold, a)
                for bridge in sorted(bridges[this.scaffold][a],key=lambda x:x.gap):
                    if not neighbour:
                        neighbour = bridge
                    print(bridge)
    return neighbour

def extend_ok(i, direction, pools):
    pool_scaffolds = {}
    cms = {}
    if direction == 1:
        print('Testing pool against next')
        print(pools[i])
        print(pools[i+direction])
    else:
        print('Testing pool against previous')
        print(pools[i+direction])
        print(pools[i])
    get_pool_scaffolds(pool_scaffolds, pools[i], cms, fillcms=True)
    get_pool_scaffolds(pool_scaffolds, pools[i+direction], cms)
    
    print(pool_scaffolds)
    genome_overlaps = get_genome_overlaps(pool_scaffolds, overlaps)
    pacbio_bridges = get_pacbio_overlaps(pool_scaffolds, overlaps)
    
    for a in pacbio_bridges:
        for b in pacbio_bridges[a]:
            print(a, b)
            for bridge in pacbio_bridges[a][b]:
                print(bridge)
    return i
    
    neighbour = get_neighbour(pools[i].rafts, pacbio_bridges, direction)
    if not neighbour:
        return i

    print("Neighbour:")
    print(neighbour)
    (this,) = pools[i].rafts
    to_remove = set()
    rafts = [raft for raft in pools[i+direction].rafts]
    for raft in rafts:
        if raft.scaffold == neighbour.hit2['scaffold']:
            print('Neighbour:', raft)
            for (scaffold, start, direction) in raft.logs:
                if direction == 1:
                    this.append(scaffold, start, direction)
                else:
                    this.prepend(scaffold, start, direction)
            print('i:{} dir:{}'.format(i, direction))
            print("Pools are currently:")
            print(pools[i+direction].rafts)
            print("And raft is:")
            print(raft)
            pools[i+direction].rafts.remove(bl)

#    for raft in to_remove:
#        pools[i+direction].rafts.remove(bl)

    print("Pools now:")
    print(pools[i])
    print(pools[i+direction])
    print("=======")

    return i

def extend_ok_pools(chromosome, overlaps):
    i = 0
    while i < (len(chromosome.pools)):
        thispool = chromosome.pools[i]
        if thispool.pooltype != 'ok':
            i += 1
            continue
        
        if i > 0:
            i = extend_ok(i, -1, chromosome.pools)

        if i < len(chromosome.pools)-1:
            i = extend_ok(i, 1, chromosome.pools)
        
        i += 1


def join_ok_pools(chromosome):
    for i in range(0, len(chromosome.pools)-1):
        pool_i = chromosome.pools[i]
        j=i+1
        while j<len(chromosome.pools):
            pool_j = chromosome.pools[j]

            if pool_i.pooltype == 'ok' and pool_j.pooltype == 'ok':
                if len(pool_i.rafts) != 1 or len(pool_j.rafts) != 1:
                    print("OK pool has more than one raft!")
                    sys.exit()
                (pool_i_raft,)=pool_i.rafts
                (pool_j_raft,)=pool_j.rafts

                extend(pool_i_raft, -1)
                extend(pool_j_raft, 0)
                for (scaffold, start, direction) in pool_j_raft.logs:
                    pool_i_raft.append(scaffold, start, direction)
                pool_j.rafts.pop()
                pool_j.pooltype=''
            else:
                break

            j += 1

    newpools = []
    for pool in chromosome:
        if pool.rafts:
            newpools.append(pool)
    
    chromosome.pools = newpoolss



def get_genome_stats(chromosomes):
    genome_stats = Stats("Genome")

    for chromosome in chromosomes:
        genome_stats += chromosomes[chromosome].stats

    print(genome_stats)

def load_map(genome):

    print("Loading map...")

    mapped_blocks = 0
    mapped_blocks_length = 0
    placed_blocks = 0
    placed_blocks_length = 0

    genome.db.execute("select distinct chromosome from scaffold_map order by chromosome")
    chromosomes = {chromosome_name: Chromosome(chromosome_name, genome) for chromosome_name, in genome.db.fetchall()}

    for name in chromosomes:
        mapped_blocks += chromosomes[name].mapped_blocks
        mapped_blocks_length += chromosomes[name].mapped_blocks_length
        placed_blocks += chromosomes[name].placed_blocks
        placed_blocks_length += chromosomes[name].placed_blocks_length
        
    print("Map blocks {} length {} of which {} blocks placed, length {}".format(mapped_blocks, mapped_blocks_length, placed_blocks, placed_blocks_length))

    return chromosomes


def reassemble(chromosomes, genome, args):

    print("Reassembly...")

    assembly = []

    for chromosome in chromosomes:
        chromosomes[chromosome].assemble(genome, args.threads)
        assembly += chromosomes[chromosome].get_scaffolds()

    get_genome_stats(chromosomes)
    
    deleted_blocks = 0
    deleted_length = 0
    for blocklist in assembly:
        for block in blocklist:
            deleted_blocks += 1
            deleted_length += block.length
            del genome.blocks[block.scaffold][block.start]

    print("Deleted {} blocks, length {}".format(deleted_blocks, deleted_length))

    left_blocks = 0
    left_length = 0
    for scaffold in genome.blocks:
        for start in genome.blocks[scaffold]:
            left_blocks += 1
            left_length += genome.blocks[scaffold][start].length
            assembly.append([genome.blocks[scaffold][start]])

    print("Left over {} blocks, length {}".format(left_blocks, left_length))

    return assembly


def output(assembly):
    
    print("Output...")
    
    genome_length = 0
    scaffolds = 0
    scaffold_lengths = []
    blocks = 0
    for blockgroup in assembly:
        scaffold_length = 0
        for block in blockgroup:
            blocks += 1
            scaffold_length += block.length
        scaffolds += 1
        genome_length += scaffold_length
        scaffold_lengths.append(scaffold_length)
    print("Blocks:{}\t".format(blocks), genome_stats(scaffold_lengths))


if __name__ == '__main__':
    
    args = get_args()

    global genome
    genome = GenomeData(args)

    chromosomes = load_map(genome)

    assembly = reassemble(chromosomes, genome, args)

    output(assembly)