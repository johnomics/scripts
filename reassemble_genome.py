#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
from statistics import mean
from os.path import isfile
import sqlite3 as db
from multiprocessing.dummy import Pool as ThreadPool

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
from Bio.Alphabet import generic_dna
from BCBio import GFF

class Chromosome:
    def __init__(self, name, cm, prev_cm=-1, next_cm=-1):
        self.name = str(name)
        self.blocks = []
        self.markers = {}
        self.add_marker(cm, prev_cm, next_cm)

    def add_marker(self, cm, prev_cm=-1, next_cm=-1):
        self.markers[cm] = Marker(cm, prev_cm, next_cm)
        
    def update_marker(self, cm, prev_cm = -1, next_cm = -1):
        if cm not in self.markers:
            self.add_marker(cm, prev_cm, next_cm)

        self.markers[cm].update_previous(prev_cm)
        self.markers[cm].update_next(next_cm)

    def add_block(self, cm, block):
        if cm not in self.cMDict:
            self.cMDict[cm] = [block]
        else:
            self.cMDict[cm].append(block)

    def __repr__(self):
        return('\n'.join([self.name, repr(self.blocks), repr(self.markers)]))

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


class BlockSet:
    def __init__(self, grouptype=''):
        self.groups = set()
        self.grouptype = grouptype

    def add(self, blocklist):
        self.groups.add(blocklist)

    def __repr__(self):
        output = ''
        for b in self.groups:
            output += repr(b) + '\n'
            output += '{}\t{}\n'.format(self.grouptype, b.length)
            output += '-----\n'

        output += '====='
        
        return output


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
            seq = sequences[self.scaffold][self.start-1:self.end]
            seq.id = '_'.join([self.scaffold, str(self.start), str(self.end)])
        else:
            seq = sequences[self.scaffold][self.end-1:self.start][::-1]
            seq.id = '_'.join([self.scaffold, str(self.start), str(self.end)])
        seq.name = seq.id
        seq.description = seq.id
        return seq

        
    def __repr__(self):
        return '{}:{}-{} ({}, {} bp)'.format(self.scaffold, str(self.start), str(self.end), str(self.cm), str(self.length))


class BlockList:
    def __init__(self, bl_id, scaffold, start):
        self.id = bl_id
        self.blocklist = []
        self.summarylist = []
        self.append(scaffold, start, 1)

    @property
    def name(self):
        summarydict = {}
        scaffolds = []
        for sb in self.summarylist:
            if sb.scaffold not in summarydict:
                summarydict[sb.scaffold] = []
                scaffolds.append(sb.scaffold)
            summarydict[sb.scaffold].append((sb.start,sb.end))
        names = []
        for scaffold in scaffolds:
            names.append('{}_{}_{}'.format(scaffold, summarydict[scaffold][0][0], summarydict[scaffold][-1][-1]))
        name = '-'.join(names)
        return name

    @property
    def start(self):
        return self.summarylist[0].start
    
    @property
    def end(self):
        return self.summarylist[-1].end

    @property
    def length(self):
        length = 0
        for sb in self.summarylist:
            length += sb.length
        return length
    
    @property
    def sequence(self):
        return sum([sb.sequence for sb in self.summarylist], Seq("", generic_dna))
    
    def append(self, scaffold, start, direction=1):
        self.blocklist = self.blocklist + [(scaffold, start, direction)]
        self.summarylist = self.summarylist + [SummaryBlock(scaffold, start, blocks[scaffold][start], direction)]
        self.update()

    def prepend(self, scaffold, start, direction=1):
        self.blocklist = [(scaffold, start, direction)] + self.blocklist
        self.summarylist = [SummaryBlock(scaffold, start, blocks[scaffold][start], direction)] + self.summarylist
        self.update()

    def update(self):
        self.collapse()
        self.chain_cms()

    def __hash__(self):
        return self.id

    def __eq__(self, other):
        return self.id == other.id

    def __iter__(self):
        return iter(self.blocklist)

    def chain_cms(self):
        cm_chain = []
        for sb in self.summarylist:
            if sb.cm != -1:
                cm_chain.append(sb.cm)
                if len(cm_chain)>1 and cm_chain[-2] == cm_chain[-1]:
                    del cm_chain[-1]

        self.cm_chain = tuple(cm_chain)
    
    def fix_chain(self):
        if check_chain(self.cm_chain):
            return True
    
        print("To fix:", self.cm_chain, "\n", self, "\n")
        return True

    def collapse_consecutive(self):
        newsummary = []
        for i in range(0, len(self.summarylist)):
            sbi = self.summarylist[i]
            if len(newsummary) == 0:
                newsummary.append(sbi)
            else:
                last = newsummary[-1]
                if last.scaffold==sbi.scaffold and last.end == sbi.start-1 and last.cm == sbi.cm:
                    last.end = sbi.end
                else:
                    newsummary.append(sbi)

        self.summarylist = newsummary

    def collapse_trios(self):
        if len(self.summarylist) < 3:
            return

        newsummary = []
        i = 0
        while i < len(self.summarylist):
            sbi = self.summarylist[i]
            if i < len(self.summarylist)-2:
                sbj = self.summarylist[i+1]
                sbk = self.summarylist[i+2]
                if (sbi.scaffold == sbj.scaffold and sbi.scaffold == sbk.scaffold and
                    sbi.end == sbj.start-1       and sbj.end == sbk.start-1 and
                    sbi.cm == sbk.cm             and sbj.cm == -1):
                        sbi.end = sbk.end
                        newsummary.append(sbi)
                        i += 3
                        continue

            newsummary.append(sbi)
            i += 1
        
        self.summarylist = newsummary
        

    def collapse(self):
        self.collapse_consecutive()
        self.collapse_trios()

    def order(self):
        if self.summarylist[0].cm != -1 and self.summarylist[-1].cm != -1 and self.summarylist[0].cm > self.summarylist[-1].cm:
            self.blocklist.reverse()
            for i in range(0, len(self.blocklist)):
                (scaffold, start, direction) = self.blocklist[i]
                self.blocklist[i] = (scaffold, start, -1)

            self.summarylist.reverse()
            for sb in self.summarylist:
                sb.start, sb.end = sb.end, sb.start
            self.chain_cms()


    def __repr__(self):
        return '\n'.join([repr(s) for s in self.summarylist])

def check_chain(chain):
    
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

def load_genome(fasta, overlaps):
    global sequences
    try:
        with open(fasta, 'r') as f:
            sequences = SeqIO.to_dict(SeqIO.parse(f, "fasta"))

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
        
    except IOError:
        print("Can't load genome from file {}!".format(fasta))
        sys.exit()


def load_annotation(gff):
    try:
        with open(gff, 'r') as g:
            for feature in g:
                if 'gene' in feature:
                    scaffold, source, gfftype, start, end, score, strand, phase, attributes = feature.rstrip().split('\t')
                    if scaffold not in sequences:
                        continue
                    featurestrand = None
                    if strand == '+':
                        featurestrand = 1
                    elif strand == '-':
                        featurestrand = -1
                    location = FeatureLocation(ExactPosition(start), ExactPosition(end), strand=featurestrand)
                    feature = SeqFeature(location, type=gfftype)
                    sequences[scaffold].features.append(feature)
    except IOError:
        print("Can't load annotation from file {}!".format(gff))
        sys.exit()

def add_overlap(rscaffold, rstart, rend, qscaffold, pcid, hittype):
    if rstart < rend:
        location = FeatureLocation(ExactPosition(rstart), ExactPosition(rend), strand=1)
    else:
        location = FeatureLocation(ExactPosition(rend), ExactPosition(rstart), strand=1)
    feature = SeqFeature(location, type="Overlap")
    sequences[rscaffold].features.append(feature)


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


def open_database(dbfile):
    try:
        if isfile(dbfile):
            conn = db.connect(dbfile)
        else:
            raise IOError
    except:
        print("Can't open database {}!".format(dbfile))
        sys.exit()
    
    return conn


def load_blocks(c, overlaps):
    blocks = {}
    block_num = 0
    genome_length = 0

    # Load map blocks from database
    for scaffold, start, end in c.execute("select scaffold, start, end from mapblocks"):
        if scaffold not in sequences:
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


def load_markers(c):
    
    chromosomes = {}
    prev_chromosome = -1
    prev_cm = -1

    for chromosome, cm in c.execute("select distinct chromosome, cm from scaffold_map order by chromosome, cm"):
        if chromosome not in chromosomes:
            chromosomes[chromosome] = Chromosome(chromosome, cm)

        if cm == -1:
            chromosomes[chromosome].add_marker(-1, -1, -1)
            continue

        if prev_chromosome != chromosome:
            prev_cm = -1
            prev_chromosome = chromosome
        
        chromosomes[chromosome].add_marker(cm, prev_cm=prev_cm)
        if prev_cm != -1:
            chromosomes[chromosome].update_marker(prev_cm, next_cm = cm)
        
        prev_cm = cm

    return chromosomes

def load_errors(errorfile):
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

def load_map(c, blocks, errors):

    mapped_blocks = 0
    mapped_blocks_length = 0
    placed_blocks = 0
    placed_blocks_length = 0
    chromosomes = load_markers(c)

    for (chromosome, chr_map) in chromosomes.items():
        for cm in sorted(chr_map.markers.keys()) + [-1]:
            cm_blocks = BlockSet()
            cm_block_id = 1
            for scaffold, start, end, length in c.execute("select scaffold, start, end, length from scaffold_map where chromosome={} and cm={} order by scaffold, start, end".format(chromosome, cm)):
                if scaffold not in sequences:
                    continue
                if (scaffold, start) in errors:
                    continue
                mapped_blocks += 1
                mapped_blocks_length += length
                blocks[scaffold][start].add_marker(chromosome,cm)

                if cm != -1:
                    placed_blocks += 1
                    placed_blocks_length += length
                cm_blocks.add(BlockList(cm_block_id, scaffold, start))
                cm_block_id += 1
            if cm != -1:
                chr_map.blocks.append(cm_blocks)

    print("Map blocks {} length {} of which {} blocks placed, length {}".format(mapped_blocks, mapped_blocks_length, placed_blocks, placed_blocks_length))

    return chromosomes

def merge_blocks(b, scaffold, block_i, chromosome, direction, markers, blocks, blocks_to_delete):
    
    block = blocks[scaffold][block_i]
    blocks_to_merge = []
    this_cm = block.cm
    
    while True:
        if direction == 1:
            block_i = blocks[scaffold][block_i].next_block
        else:
            block_i = blocks[scaffold][block_i].prev_block

        if block_i == 0:
            break
        
        block_i_chr = blocks[scaffold][block_i].chromosome
        block_i_cm = blocks[scaffold][block_i].cm
        
        if block_i_chr == '0' or block_i_cm == -1 and block_i_chr == chromosome:
            blocks_to_merge.append(block_i)
            continue

        if block_i_chr != chromosome:
            break

        merge = False
        
        if block_i_cm == this_cm:
            merge = True
        
        if block_i_cm == markers[this_cm].next_cm:
            this_cm = markers[this_cm].next_cm
            merge = True
        
        if block_i_cm == markers[this_cm].prev_cm:
            this_cm = markers[this_cm].prev_cm
            merge = True

        if merge:
            blocks_to_merge.append(block_i)
            for merge_block in blocks_to_merge:
                if scaffold not in blocks_to_delete:
                    blocks_to_delete[scaffold] = {}
                blocks_to_delete[scaffold][merge_block] = True
                
                if direction == 1:
                    b.append(scaffold, merge_block)
                else:
                    b.prepend(scaffold, merge_block)
            
            blocks_to_merge=[]
            continue

        break


def collapse_neighbours(chromosome, blocks):

    blocks_to_delete = {}
    
    for blockset in chromosome.blocks:
        blocklists_to_delete = []
        for b in blockset.groups:
            (scaffold, start, direction) = b.blocklist[0]

            if scaffold in blocks_to_delete and start in blocks_to_delete[scaffold]:
                blocklists_to_delete.append(b)
                continue
            
            merge_blocks(b, scaffold, start, chromosome.name, 1, chromosome.markers, blocks, blocks_to_delete)
            merge_blocks(b, scaffold, start, chromosome.name, -1, chromosome.markers, blocks, blocks_to_delete)
            
        for b in blocklists_to_delete:
            blockset.groups.remove(b)

    newblocks = []
    for blockset in chromosome.blocks:
        if blockset.groups:
            newblocks.append(blockset)

    chromosome.blocks = newblocks


def refine_groups(chromosome):
    """If a well-ordered block exists in a group,
       split it off into its own group
    """
    newblocks = []
    for blockset in chromosome.blocks:
        newblocksets = [BlockSet()]
        for b in blockset.groups:
            b.order()
            if b.summarylist[0].cm != b.summarylist[-1].cm:
                newblocksets.append(BlockSet())
                newblocksets[-1].groups.add(b)
            else:
                newblocksets[0].groups.add(b)

        if not newblocksets[0].groups:
            del newblocksets[0]

        newblocks += newblocksets

    chromosome.blocks = newblocks

def extend_groups(chromosome):
    for blockset in chromosome.blocks:
        for b in blockset.groups:
            extend_to_ends(b)


class GroupStat:
    def __init__(self, count=0, length=0):
        self.count = count
        self.length = length

def classify_groups(chromosome):

    cm_chains = {}

    for blockset in chromosome.blocks:
        for b in blockset.groups:
            b.fix_chain()
            
            if b.cm_chain not in cm_chains:
                cm_chains[b.cm_chain] = {}
                cm_chains[b.cm_chain]['count'] = 0
                cm_chains[b.cm_chain]['type'] = ''
            
            cm_chains[b.cm_chain]['count'] += 1
    
    group_stats = {'orient':GroupStat(), 'order':GroupStat(), 'overlap':GroupStat(), 'error':GroupStat(), 'ok':GroupStat()}
    
    for cms in sorted(cm_chains):
        
        group_type = ''
        if len(cms) == 1:
            if cm_chains[cms]['count'] == 1:
                group_type = 'orient'
            else:
                group_type = 'order'
        else:
            if cm_chains[cms]['count'] > 1:
                group_type = 'overlap'
            else:
                if check_chain(cms):
                    group_type = 'ok'
                else:
                    group_type = 'error'
        
        if group_type != '':
            group_stats[group_type].count += 1
            cm_chains[cms]['type'] = group_type

    block_num = 0
    for blockset in chromosome.blocks:
        for b in blockset.groups:
            block_num += len(b.blocklist)
            gt = cm_chains[b.cm_chain]['type']
            blockset.grouptype = gt
            group_stats[gt].length += b.length

    print("{:3s}\t{:6d}".format(chromosome.name, block_num),end='')
    for gt in 'ok', 'orient', 'order', 'overlap', 'error':
        print("\t{:6d}\t{:10d}".format(group_stats[gt].count, group_stats[gt].length), end='')
    print()


def identify_groups(chromosome, blocks):
    collapse_neighbours(chromosome, blocks)
    refine_groups(chromosome)
    extend_groups(chromosome)
    classify_groups(chromosome)


def call_finisher(dirname):
    with open(dirname + '/finisherSC.log', 'w') as logfile:
        subprocess.call(['python', '/whale-data/jd626/bin/finishingTool/finisherSC.py', '-o', 'contigs.fasta_improved3.fasta', '.', '/whale-data/jd626/bin/MUMmer/elephant/MUMmer3.23'], stdout=logfile, stderr=logfile, cwd=dirname)

def build_groups(chromosome, threads):
    dirnames = []
    for block in chromosome.blocks:
        for b in block.groups:
            for sb in b.summarylist:
                if sb.cm == -1:
                    continue
                chrstring = 'chr' + str(chromosome.name)
                dirname = 'reassembly/' + chrstring + '/' + chrstring + '_' + str(sb.cm)
                filename = dirname + '/contigs.fasta'
                linkname = dirname + '/raw_reads.fasta'
                if not os.path.exists(dirname):
                    os.makedirs(dirname, exist_ok=True)
                    os.symlink("/whale-data/jd626/hmel_pacbio/LoRDEC/heliconius_melpomene_melpomene.subreads.LoRDEC.with_tail.fasta", linkname)
                    dirnames.append(dirname)
                with open(filename, 'a') as f:
                    SeqIO.write(sb.sequence, f, "fasta")

    pool = ThreadPool(threads)
    pool.map(call_finisher, dirnames)
    pool.close()
    pool.join()


def extend(b, item):
    
    first_scaffold, first_start, direction = b.blocklist[item]
    ext_start = first_start
    extend = 0
    starts_to_extend = []
    chromosome = blocks[first_scaffold][first_start].chromosome

    while 1:
        if item == 0 and direction == 1 or item == -1 and direction == -1:
            ext_start = blocks[first_scaffold][ext_start].prev_block
        else:
            ext_start = blocks[first_scaffold][ext_start].next_block

        if ext_start == 0:
            break

        if blocks[first_scaffold][ext_start].cm != -1 or blocks[first_scaffold][ext_start].chromosome != 0 and blocks[first_scaffold][ext_start].chromosome != chromosome:
            break

        starts_to_extend.append(ext_start)

    if starts_to_extend:
        if item == 0:
            for start in starts_to_extend:
                b.prepend(first_scaffold, start, direction)
        else:
            for start in starts_to_extend:
                b.append(first_scaffold, start, direction)


def extend_to_ends(b):
    extend(b, 0)
    extend(b, -1)
    

def build_scaffolds(chromosome, blocks):

    scaffolds = []
    groups = {}

    for blockset in chromosome.blocks:
        setnum = len(blockset.groups)
        if setnum not in groups:
            groups[setnum] = 0
        groups[setnum] += 1

        for b in blockset.groups:
            scaffoldlist = []
            
            extend_to_ends(b)

            for scaffold, start, direction in b.blocklist:
                scaffoldlist.append(blocks[scaffold][start])
            scaffolds.append(scaffoldlist)

    return scaffolds, groups

def get_block_scaffolds(scaffolds, blockset, cms=None, fillcms=False):
    for blocks in blockset.groups:
        for sb in blocks.summarylist:
            if fillcms and sb.cm != -1 and sb.cm not in cms:
                cms[sb.cm] = 0
            if sb.cm in cms and sb.scaffold not in scaffolds:
                scaffolds[sb.scaffold] = {'start': blocks.start, 'end': blocks.end}


def get_group_scaffolds(i, blocks):
    group_scaffolds = {}
    cms = {}
    print(blocks[i])
    get_block_scaffolds(group_scaffolds, blocks[i], cms, fillcms=True)
    if i > 0:
        get_block_scaffolds(group_scaffolds, blocks[i-1], cms)
    if i < len(blocks)-1:
        get_block_scaffolds(group_scaffolds, blocks[i+1], cms)

    return group_scaffolds, cms

def bridge_with_genome_overlaps(a, b, overlaps):
    found = False
    statement = "select * from genome_overlaps where rscaffold=\"{}\" and qscaffold=\"{}\"".format(a, b)
    for rstart, rend, qstart, qend, rhitlen, qhitlen, pcid, rseqlen, qseqlen, rpccov, qpccov, rscaffold, qscaffold, hittype in overlaps.execute(statement):
        found = True
        print(a, len(sequences[a]), rstart, rend, rhitlen, rpccov, b, len(sequences[b]), qstart, qend, qhitlen, qpccov, pcid, hittype, sep='\t')
    
    return found


def print_genome_overlaps(scaffolds, overlaps):
    found = False
    for a in scaffolds:
        for b in scaffolds:
            if a == b:
                continue
            abfound = bridge_with_genome_overlaps(a, b, overlaps)
            if abfound:
                found = True
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
        bstart, bend = order_hit_ends(scaffolds[qscaffold]['start'], scaffolds[qscaffold]['end'])
        if (gstart > bstart and gend < bend
           or gstart > bend or gend < bstart):
            continue
        if rscaffold not in pacbio_scaffolds:
            pacbio_scaffolds[rscaffold] = {}
            pacbio_lengths[rscaffold] = rseqlen
        if qscaffold not in pacbio_scaffolds[rscaffold]:
            pacbio_scaffolds[rscaffold][qscaffold] = []
        pacbio_scaffolds[rscaffold][qscaffold].append((qstart, qend, qhitlen, rstart, rend, pcid))

    return pacbio_scaffolds, pacbio_lengths


def get_pacbio_bridges(scaffolds, pacbio_scaffolds, pacbio_lengths):
    bridges = {}
    for pbs in pacbio_scaffolds:
        if len(pacbio_scaffolds[pbs]) <= 1:
            continue
        pacbio_order = {}
        pairs = {tuple(sorted([a,b])) for a in pacbio_scaffolds[pbs] for b in pacbio_scaffolds[pbs] if a != b}
        for gs in pacbio_scaffolds[pbs]:
            for (qstart, qend, qhitlen, rstart, rend, pcid) in sorted(pacbio_scaffolds[pbs][gs], key=lambda x:x[2], reverse=True):
                if rstart > rend:
                    qstart, qend = qend, qstart
                    rstart, rend = rend, rstart
                ignore = False
                for start, hit in pacbio_order.items():
                    if gs == hit['scaffold'] and qstart >= hit['gstart'] and qend <= hit['gend']:
                        ignore = True
                if not ignore:
                    pacbio_order[rstart] = {'scaffold': gs, 'pbstart': rstart, 'pbend': rend, 'gstart': qstart, 'gend': qend, 'pcid':pcid}

        for s1, hit1 in sorted(pacbio_order.items()):
            for s2, hit2 in sorted(pacbio_order.items()):
                if s2 <= s1:
                    continue
                print('''{0} ({1:6d} bp) {2} - {3} joins
        {4:16s} ({5:7d} bp) {6:7d} - {7:7d} ({2:7d} - {8:7d}; {9:5.2f}%) and
        {10:16s} ({11:7d} bp) {12:7d} - {13:7d} ({14:7d} - {3:7d}; {15:5.2f}%)'''.format(
                    pbs, pacbio_lengths[pbs], hit1['pbstart'], hit2['pbend'],
                    hit1['scaffold'], len(sequences[hit1['scaffold']]),
                    hit1['gstart'], hit1['gend'], hit1['pbend'], hit1['pcid'],
                    hit2['scaffold'], len(sequences[hit2['scaffold']]),
                    hit2['gstart'], hit2['gend'], hit2['pbstart'], hit2['pcid']
                ), end='')
                if hit1['pbend'] < hit2['pbstart']:
                    print("\tFill gap with {} {}-{}".format(pbs, hit1['pbend']+1, hit2['pbstart']-1))
                elif hit1['pbstart'] < hit2['pbstart'] and hit1['pbend'] > hit2['pbend']:
                    print("\tSecond hit contained in first hit")
                elif hit2['pbstart'] < hit1['pbstart'] and hit2['pbend'] > hit1['pbend']:
                    print("\tFirst hit contained in second hit")
                else:
                    print("\tOverlap of {} bases".format(hit1['pbend'] - hit2['pbstart'] + 1))

    return bridges


def print_pacbio_overlaps(scaffolds, overlaps):
    pacbio_scaffolds, pacbio_lengths = get_pacbio_hits(scaffolds, overlaps)
    
    if not pacbio_scaffolds:
        print("No PacBio overlaps found")
        return

    bridges = get_pacbio_bridges(scaffolds, pacbio_scaffolds, pacbio_lengths)


def view_groups(chromosome, overlaps):
    for i in range(len(chromosome.blocks)):
        group_scaffolds, cms = get_group_scaffolds(i, chromosome.blocks)

        if len(cms) > 1:
            continue

#        print(sorted(cms.keys()))
#        print([(g, len(sequences[g])) for g in group_scaffolds])

#        print_genome_overlaps(group_scaffolds, overlaps)
        print_pacbio_overlaps(group_scaffolds, overlaps)
        
        print("=======")


def get_hits(sl, group):
    for blockset in group.groups:
        for sb in blockset.summarylist:
            statement = "select * from genome_overlaps where rscaffold=\"{}\" and qscaffold=\"{}\"".format(sl.scaffold, sb.scaffold)
#            statement = '''select * from overlaps
#                           where rscaffold=\"{}\" and (rstart <= {} and rend >= {} or rstart <= {} and rend >= {}) and 
#                           qscaffold=\"{}\" and (rstart <= {} and rend >= {} or rstart <= {} and rend >= {})'''.format(
#                               sl.scaffold, sl.start, sl.start, sl.end, sl.end,
#                               sb.scaffold, sb.start, sb.start, sb.end, sb.end)
            for rstart, rend, qstart, qend, rhitlen, qhitlen, pcid, rpccov, qpccov, rscaffold, qscaffold, hittype in overlaps.execute(statement):
                rout = "{} (length {}) {}-{} ({} bp, {} covered)".format(sl.scaffold, len(sequences[sl.scaffold]), rstart, rend, rhitlen, rpccov)
                qout = "{} (length {}) {}-{} ({} bp, {} covered)".format(sb.scaffold, len(sequences[sb.scaffold]), qstart, qend, qhitlen, qpccov)
                print(rout, qout, pcid, hittype, sep='\t')
    

def extend_ok_groups(chromosome, overlaps):
    
    for i in range(0, len(chromosome.blocks)):
        thisgroup = chromosome.blocks[i]
        if thisgroup.grouptype != 'ok':
            continue

        (thisblockset,) = thisgroup.groups
        prevsl = thisblockset.summarylist[0]
        nextsl = thisblockset.summarylist[-1]
        
        print(prevsl.cm, nextsl.cm)

        if i > 0:
            prevgroup = chromosome.blocks[i-1]
            print(prevgroup)
            get_hits(prevsl, prevgroup)

        print(thisgroup)
            
        if i < len(chromosome.blocks)-1:
            nextgroup = chromosome.blocks[i+1]
            print(nextgroup)
            get_hits(nextsl, nextgroup)

    newblocks = []
    for blockset in chromosome.blocks:
        if blockset.groups:
            newblocks.append(blockset)
    
    chromosome.blocks = newblocks
    

def join_ok_groups(chromosome):
    for i in range(0, len(chromosome.blocks)-1):
        bi = chromosome.blocks[i]
        j=i+1
        while j<len(chromosome.blocks):
            bj = chromosome.blocks[j]

            if bi.grouptype == 'ok' and bj.grouptype == 'ok':
                if len(bi.groups) != 1 or len(bj.groups) != 1:
                    print("OK group has more than one blocklist!")
                    sys.exit()
                (bi_bl,)=bi.groups
                (bj_bl,)=bj.groups

                extend(bi_bl, -1)
                extend(bj_bl, 0)
                for (scaffold, start, direction) in bj_bl.blocklist:
                    bi_bl.append(scaffold, start, direction)
                bj.groups.pop()
                bj.grouptype=''
            else:
                break

            j += 1

    newblocks = []
    for blockset in chromosome.blocks:
        if blockset.groups:
            newblocks.append(blockset)
    
    chromosome.blocks = newblocks


def assemble_chromosome(chromosome, blocks, threads, overlaps):

    identify_groups(chromosome, blocks)

    view_groups(chromosome, overlaps)

    scaffolds, groups = build_scaffolds(chromosome, blocks)

    return scaffolds, groups


def reassemble(chromosomes, blocks, threads, overlaps):
    genome = []

    groups = {}
    print("Chromosome\tBlocks\tOK\tOrient\tOrder\tOverlap\tError")

    for chromosome in chromosomes:
        chromosome_blocks, chromosome_groups = assemble_chromosome(chromosomes[chromosome], blocks, threads, overlaps)
        genome += chromosome_blocks
        
        for i in chromosome_groups:
            if i not in groups:
                groups[i] = 0
            groups[i] += chromosome_groups[i]
    
    for i in groups:
        print(i, groups[i])
    
    
    deleted_blocks = 0
    deleted_length = 0
    for blocklist in genome:
        for block in blocklist:
            deleted_blocks += 1
            deleted_length += block.length
            del blocks[block.scaffold][block.start]

    print("Deleted {} blocks, length {}".format(deleted_blocks, deleted_length))

    left_blocks = 0
    left_length = 0
    for scaffold in blocks:
        for start in blocks[scaffold]:
            left_blocks += 1
            left_length += blocks[scaffold][start].length
            genome.append([blocks[scaffold][start]])

    print("Left over {} blocks, length {}".format(left_blocks, left_length))

    return genome


def output(genome):
    
    genome_length = 0
    scaffolds = 0
    scaffold_lengths = []
    blocks = 0
    for blockgroup in genome:
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


    print("Loading genome...")
    overlap_db = open_database(args.overlaps)
    overlaps = overlap_db.cursor()

    load_genome(args.fasta, overlaps)
    
    print("Loading annotation...")
    load_annotation(args.gff)
    
    
    conn = open_database(args.database)
    c = conn.cursor()

    print("Loading blocks...")
    blocks = load_blocks(c, overlaps)

    print("Loading errors...")
    errors = load_errors(args.errors)

    print("Loading map...")

    chromosomes = load_map(c, blocks, errors)

    print("Reassembly...")
    genome = reassemble(chromosomes, blocks, args.threads, overlaps)

    print("Output...")
    output(genome)
