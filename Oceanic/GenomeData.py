#!/usr/bin/env python3

import sys
import gzip
import re
from os.path import isfile
import sqlite3 as sql
from collections import defaultdict

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
from BCBio import GFF

from . import Stats

class GenomeData:
    def __init__(self, args):
        self.haplotypes = {}

        print("Loading genome...")
        
        self.sequences = self.load_genome(args.fasta)
        
        print("Loading annotation...")
        self.load_annotation(args.gff)
        
        conn = self.open_database(args.database)
        self.db = conn.cursor()
        
        print("Loading blocks...")
        self.blocks = self.load_blocks()
        
        print("Loading errors...")
        self.errors = self.load_errors(args.errors)
        
        self.revised, self.revised_fasta, self.revised_db, self.revised_conn = self.open_revised(args.revised)

        self.gapnum = 1
        
        if args.haplomerger:
            print("Loading haplotypes...")
            self.haplotypes = self.load_haplotypes(args.haplomerger)
        else:
            self.haplotypes = {}
        
    def open_revised(self, revised):
        fasta = None
        if revised:
            fasta = open(revised+".fa", 'w')
            rconn = sql.connect(revised+".db")
            db = rconn.cursor()
            db.execute('drop table if exists scaffold_map')
            db.execute('''create table scaffold_map
                         (chromosome integer, cm real, scaffold text, start integer, end integer, length integer)''')

        return revised, fasta, db, rconn

    def open_database(self, dbfile):
        try:
            if isfile(dbfile):
                conn = sql.connect(dbfile)
            else:
                raise IOError
        except IOError:
            print("Can't open database {}!".format(dbfile))
            sys.exit()
    
        return conn

    def load_genome(self, fasta):
        try:
            with open(fasta, 'r') as f:
                sequences = SeqIO.to_dict(SeqIO.parse(f, "fasta"))
        except IOError:
            print("Can't load genome from file {}!".format(fasta))
            sys.exit()

        print("Original stats:", Stats.genome([len(sequences[scaffold]) for scaffold in sequences]))

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

    def load_blocks(self):
        blocks = defaultdict(lambda:defaultdict(Block))
        block_num = 0
        genome_length = 0

        # Load map blocks from database
        for scaffold, start, end in self.db.execute("select scaffold, start, end from scaffold_map"):
            if scaffold not in self.sequences:
                continue
            if end < start:
                continue
            block_num += 1
            genome_length += end - start + 1
            blocks[scaffold][start] = Block(scaffold, start, end)
        
        for scaffold in self.sequences:
            if scaffold not in blocks:
                end = len(self.sequences[scaffold])
                block_num += 1
                genome_length += end
                blocks[scaffold][1] = Block(scaffold, 1, end)

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
        except IOError:
            print("Can't open errors file!")
            sys.exit()
    
        return errors

    def add_gap(self, length = 100):
        newgapname = 'Gap' + str(self.gapnum)
        self.blocks[newgapname][1] = Block(newgapname, 1, length)
        self.gapnum += 1
        return self.blocks[newgapname][1]
    
    def load_haplotypes(self, haplomerger):
        haplotypes = {}
        self.load_haplotype(haplotypes, haplomerger, "A")
        self.load_haplotype(haplotypes, haplomerger, "B")
        return haplotypes
    
    def load_haplotype(self, haplotypes, haplomerger, hap):
        hapfile = haplomerger + "/assembly_hap" + hap + ".fa.gz"
        try:
            if isfile(hapfile):
                with gzip.open(hapfile, 'rt') as hf:
                    for line in hf:
                        if line.startswith(">"):
                            f = line[1:].split('_')
                            name = f[0] + '_' + f[1] + '_' + f[2]
                            haplotypes[name] = hap
        except IOError:
            print("Can't open haplotype file!", hapfile)
            sys.exit()
    

class Block:
    def __init__(self, scaffold, start, end, prev_block=0, next_block=0, chromosome=0, cm=-1):
        self.scaffold = scaffold
        self.start = start
        self.end = end
        self.prev_block = prev_block
        self.next_block = next_block
        self.chromosome = str(chromosome)
        self.cm = cm
        self.contained = None
    
    @property
    def length(self):
        return self.end - self.start + 1
    
    def __repr__(self):
        return "{}:{}-{} ({}, {}) [{}, {}]".format(self.scaffold, self.start, self.end, self.prev_block, self.next_block, self.chromosome, self.cm)

    def add_marker(self, chromosome, cm):
        self.chromosome = str(chromosome)
        self.cm = cm

if __name__ == '__main__':
    pass