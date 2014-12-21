#!/usr/bin/env python3

import sys
import argparse
from os.path import isfile
import sqlite3 as db
from Bio import SeqIO

class Chromosome:
    def __init__(self, name, cm, prev_cm=-1, next_cm=-1):
        self.name = name
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
        self.length = end - start + 1
        self.prev_block = prev_block
        self.next_block = next_block
        self.chromosome = chromosome
        self.cm = cm
    
    def __repr__(self):
        return "{}:{}-{} ({}, {})".format(self.scaffold, self.start, self.end, self.prev_block, self.next_block)

    def add_marker(self, chromosome, cm):
        self.chromosome = chromosome
        self.cm = cm


class SummaryBlock:
    def __init__(self, scaffold, start, block):
        self.scaffold = scaffold
        self.start = start
        self.end = block.end
        self.length = max(self.start,self.end)-min(self.start,self.end)+1
        self.cm = block.cm

    def __repr__(self):
        return '{}:{}-{} ({}, {} bp)'.format(self.scaffold, str(self.start), str(self.end), str(self.cm), str(self.length))

    def update_length(self):
        self.length = max(self.start,self.end)-min(self.start,self.end)+1


class BlockList:
    def __init__(self, bl_id, scaffold, start):
        self.blocklist = [(scaffold, start)]
        self.summarylist = [SummaryBlock(scaffold, start, blocks[scaffold][start])]
        self.id = bl_id
        self.get_length()
        self.chain_cms()
    
    def append(self, scaffold, start):
        self.blocklist = self.blocklist + [(scaffold, start)]
        self.summarylist = self.summarylist + [SummaryBlock(scaffold, start, blocks[scaffold][start])]
        self.update()

    def prepend(self, scaffold, start):
        self.blocklist = [(scaffold, start)] + self.blocklist
        self.summarylist = [SummaryBlock(scaffold, start, blocks[scaffold][start])] + self.summarylist
        self.update()

    def update(self):
        self.collapse()
        self.chain_cms()
        self.get_length()

    def __hash__(self):
        return self.id

    def __eq__(self, other):
        return self.id == other.id

    def __iter__(self):
        return iter(self.blocklist)

    def get_length(self):
        self.length = 0
        for sb in self.summarylist:
            self.length += sb.length

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
                    last.update_length()
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
                        sbi.update_length()
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
        -e errors file''')

    parser.add_argument('-d', '--database', type=str, required=True)
    parser.add_argument('-f', '--fasta', type=str, required=True)
    parser.add_argument('-g', '--gff', type=str, required=True)
    parser.add_argument('-e', '--errors', type=str)

    return parser.parse_args()


def open_database(dbfile):
    try:
        if isfile(dbfile):
            conn = db.connect(dbfile)
        else:
            raise IOError
    except:
        print("Can't open database!")
        sys.exit()
    
    return conn


def load_blocks(c):
    blocks = {}
    block_num = 0
    genome_length = 0

    # Load map blocks from database
    for scaffold, start, end in c.execute("select scaffold, start, end from mapblocks"):
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

        if prev_chromosome != chromosome:
            prev_cm = -1
            prev_chromosome = chromosome
            
        if chromosome not in chromosomes:
            chromosomes[chromosome] = Chromosome(chromosome, cm)
        
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

    genome_length = 0
    block_num = 0
    chromosomes = load_markers(c)

    for (chromosome, chr_map) in chromosomes.items():
        for (cm, marker) in sorted(chr_map.markers.items()):
            cm_blocks = set()
            cm_block_id = 1
            for scaffold, start, end, length in c.execute("select scaffold, start, end, length from scaffold_map where chromosome={} and cm={} order by scaffold, start, end".format(chromosome, cm)):
                if (scaffold, start) in errors:
                    continue
                block_num += 1
                genome_length += length
                blocks[scaffold][start].add_marker(chromosome,cm)
                cm_blocks.add(BlockList(cm_block_id, scaffold, start))
                cm_block_id += 1
            chr_map.blocks.append(cm_blocks)

    print("Map blocks {} length {}".format(block_num, genome_length))

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
        if block_i_chr == 0:
            blocks_to_merge.append(block_i)
            continue

        if block_i_chr != chromosome:
            break

        block_i_cm = blocks[scaffold][block_i].cm

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
        for b in blockset:
            (scaffold, start) = b.blocklist[0]

            if scaffold in blocks_to_delete and start in blocks_to_delete[scaffold]:
                blocklists_to_delete.append(b)
                continue
            
            merge_blocks(b, scaffold, start, chromosome.name, 1, chromosome.markers, blocks, blocks_to_delete)

            merge_blocks(b, scaffold, start, chromosome.name, -1, chromosome.markers, blocks, blocks_to_delete)
            
        for b in blocklists_to_delete:
            blockset.remove(b)

    newblocks = []
    for blockset in chromosome.blocks:
        if blockset:
            newblocks.append(blockset)

    chromosome.blocks = newblocks


def refine_groups(chromosome):
    newblocks = []
    for blockset in chromosome.blocks:
        newblocksets = [set()]
        for b in blockset:
            b.order()
            if b.summarylist[0].cm != b.summarylist[-1].cm:
                newblocksets.append(set())
                newblocksets[-1].add(b)
            else:
                newblocksets[0].add(b)

        if not newblocksets[0]:
            del newblocksets[0]

        newblocks += newblocksets

    chromosome.blocks = newblocks
    

def assemble_chromosome(chromosome, blocks):
    scaffolds = []
    
    collapse_neighbours(chromosome, blocks)

    refine_groups(chromosome)

    cm_chains = {}
    groups = {}

    for blockset in chromosome.blocks:
        scaffoldlist = []
        setnum = len(blockset)
        if setnum not in groups:
            groups[setnum] = 0
        groups[setnum] += 1

        for b in blockset:
#            print(b)
#            print('-----')
            b.fix_chain()
            
            if b.cm_chain not in cm_chains:
                cm_chains[b.cm_chain] = {}
                cm_chains[b.cm_chain]['count'] = 0
                cm_chains[b.cm_chain]['type'] = ''
            
            cm_chains[b.cm_chain]['count'] += 1

            for scaffold, start in b.blocklist:
                scaffoldlist.append(blocks[scaffold][start])
#        print('=====')
        scaffolds.append(scaffoldlist)

    group_types = {'orient':0, 'order':0, 'overlap':0, 'error':0, 'ok':0}
    
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
            group_types[group_type] += 1
            cm_chains[cms]['type'] = group_type

    type_lengths = {'orient':0, 'order':0, 'overlap':0, 'error':0, 'ok':0}

    block_num = 0
    for blockset in chromosome.blocks:
        for b in blockset:
            block_num += len(b.blocklist)
            gt = cm_chains[b.cm_chain]['type']
            type_lengths[gt] += b.length
            print(b)
            print(gt, b.length)
            print('-----')

        print('=====')

    print("{:3d}\t{:6d}".format(chromosome.name, block_num),end='')
    for gt in 'ok', 'orient', 'order', 'overlap', 'error':
        print("\t{:6d}\t{:10d}".format(group_types[gt], type_lengths[gt]), end='')
    print()

    return scaffolds, groups


def reassemble(chromosomes, blocks):
    genome = []

    groups = {}
    print("Chromosome\tBlocks\tOK\tOrient\tOrder\tOverlap\tError")

    for chromosome in chromosomes:
        chromosome_blocks, chromosome_groups = assemble_chromosome(chromosomes[chromosome], blocks)
        genome += chromosome_blocks
        
        for i in chromosome_groups:
            if i not in groups:
                groups[i] = 0
            groups[i] += chromosome_groups[i]
    
#    for i in groups:
#        print(i, groups[i])
    
    
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
            left_length += blocks[scaffold][start].end - start + 1
            genome.append([blocks[scaffold][start]])

    print("Left over {} blocks, length {}".format(left_blocks, left_length))

    return genome


def output(genome):
    
    genome_length = 0
    scaffolds = 0
    blocks = 0
    
    for i, blockgroup in enumerate(genome):
        # Does the blocklist overlap any genes?
        # Construct FASTA for blocklist from original genome
        scaffolds += 1
        for block in blockgroup:
            blocks += 1
            genome_length += block.end - block.start + 1
            
    print("Scaffolds:{}\tBlocks:{}\tLength:{}".format(scaffolds, blocks, genome_length))



if __name__ == '__main__':
    
    args = get_args()
    
    conn = open_database(args.database)
    c = conn.cursor()

    blocks = load_blocks(c)
    errors = load_errors(args.errors)
    chromosomes = load_map(c, blocks, errors)
    genome = reassemble(chromosomes, blocks)
    output(genome)
