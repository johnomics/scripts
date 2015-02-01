#!/usr/bin/env python3

import sys
from os.path import isfile
import sqlite3 as sql

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
from BCBio import GFF

from . import Stats

class GenomeData:
    def __init__(self, args):
        self.haplotypes = {}
        
        print("Loading genome...")
        overlap_db = self.open_database(args.overlaps)
        self.overlaps = overlap_db.cursor()
        
        self.sequences = self.load_genome(args.fasta)
        
        print("Loading annotation...")
        self.load_annotation(args.gff)
        
        conn = self.open_database(args.database)
        self.db = conn.cursor()
        
        print("Loading blocks...")
        self.blocks = self.load_blocks()
        
        print("Loading errors...")
        self.errors = self.load_errors(args.errors)


    def open_database(self, dbfile):
        try:
            if isfile(dbfile):
                conn = sql.connect(dbfile)
            else:
                raise IOError
        except:
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

        delete_scaffolds = []
        for scaffold in sequences:
            statement = "select hittype from genome_overlaps where rscaffold=\"{}\"".format(scaffold)
            skip = 0
            for (hittype, ) in self.overlaps.execute(statement):
                if hittype == '[CONTAINED]':
                    skip = 1
            if skip:
                delete_scaffolds.append(scaffold)
        
        for scaffold in delete_scaffolds:
            del sequences[scaffold]

        print("Remove contained scaffolds:", Stats.genome([len(sequences[scaffold]) for scaffold in sequences]))

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
        blocks = {}
        block_num = 0
        genome_length = 0

        # Load map blocks from database
        for scaffold, start, end in self.db.execute("select scaffold, start, end from mapblocks"):
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

if __name__ == '__main__':
    pass