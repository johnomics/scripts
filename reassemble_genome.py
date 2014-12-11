#!/usr/bin/env python3

import sys
import argparse
from os.path import isfile
import sqlite3 as db
from Bio import SeqIO

class Chromosome:
    def __init__(self, chromosome, cm, prev_cm=-1, next_cm=-1):
        self.chromosome = chromosome
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


def load_map(c, blocks):

    genome_length = 0
    block_num = 0

    chromosomes = load_markers(c)
    
    for (chromosome, chr_map) in chromosomes.items():
        for (cm, marker) in sorted(chr_map.markers.items()):
            cm_blocks = set()
            for scaffold, start, end, length in c.execute("select scaffold, start, end, length from scaffold_map where chromosome={} and cm={} order by scaffold, start, end".format(chromosome, cm)):
                block_num += 1
                genome_length += length
                blocks[scaffold][start].add_marker(chromosome,cm)
                cm_blocks.add(((scaffold,start),))
            chr_map.blocks.append(cm_blocks)

    print("Map blocks {} length {}".format(block_num, genome_length))

    return chromosomes



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

def get_args():
    parser = argparse.ArgumentParser(description='''Output new FASTA file based on linkage map.
    
        -d database
        -f FASTA file
        -g GFF file''')

    parser.add_argument('-d', '--database', type=str, required=True)
    parser.add_argument('-f', '--fasta', type=str, required=True)
    parser.add_argument('-g', '--gff', type=str, required=True)

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




def assemble_chromosome(chromosome, blocks):
    chromosome_blocks = []
    for blockgroup in chromosome.blocks:
        for blockset in blockgroup:
            blocklist = []
            for block in blockset:
                (scaffold, start) = block
                blocklist.append(blocks[scaffold][start])
            chromosome_blocks.append(blocklist)

    return chromosome_blocks


def reassemble(chromosomes, blocks):
    genome = []

    for chromosome in chromosomes:
        chromosome_blocks = assemble_chromosome(chromosomes[chromosome], blocks)
        genome += chromosome_blocks
    
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
    chromosomes = load_map(c, blocks)
    genome = reassemble(chromosomes, blocks)

    output(genome)
