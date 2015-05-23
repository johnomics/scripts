#!/usr/bin/python3 -u

import os
import sys
import pathlib
import argparse
import sqlite3 as sql
from Bio import SeqIO

def load_genome(fasta):
    try:
        with open(fasta, 'r') as f:
            sequences = SeqIO.to_dict(SeqIO.parse(f, "fasta"))
    except IOError:
        print("Can't load genome from file {}!".format(fasta))
        sys.exit()

    return sequences
    
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


def output_chromosomes(genome, db, basename):
    
    dirname = "{}_chromosomes".format(basename)
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    db.execute('select distinct chromosome from scaffold_map')
    chromosomes = sorted([x[0] for x in db.fetchall()])
    del chromosomes[0]  # Remove chr 0
    
    chromfh = {}
    for chrom in chromosomes:
        chromfh[chrom] = open('{}/{}_chromosome_{}.fasta'.format(dirname,basename,chrom), 'w')
    
    misassembled = open('{}/{}_misassembled.fasta'.format(dirname,basename), 'w')
    unassigned = open('{}/{}_unassigned.fasta'.format(dirname,basename), 'w')
    removed = open('{}/{}_removed.txt'.format(dirname,basename), 'w')
    db.execute('select distinct scaffold from scaffold_map')
    scaffolds = sorted([x[0] for x in db.fetchall()])
    
    for scaffold in scaffolds:
        db.execute('select distinct chromosome from scaffold_map where scaffold=\"{}\"'.format(scaffold))
        scfchroms = sorted([x[0] for x in db.fetchall()])
        if scfchroms[0] == 0:
            del scfchroms[0]
        if scaffold not in genome:
            removed.write('{}\t{}\n'.format(scaffold, scfchroms))
            continue

        if not scfchroms:
            SeqIO.write(genome[scaffold], unassigned, "fasta")
        elif len(scfchroms)>1:
            SeqIO.write(genome[scaffold], misassembled, "fasta")
        else:
            SeqIO.write(genome[scaffold], chromfh[scfchroms[0]], "fasta")
    
    for chrom in chromfh:
        chromfh[chrom].close()

    
    
def get_args():
    parser = argparse.ArgumentParser(description='''Output chromosomes from genome based on linkage map
        -f fasta
        -d database
        ''')

    parser.add_argument('-f', '--fasta', type=str, required=True)
    parser.add_argument('-d', '--database', type=str, required=True)
    return parser.parse_args()

if __name__ == '__main__':
    
    args = get_args()

    basename = os.path.splitext(os.path.basename(args.fasta))[0]

    genome = load_genome(args.fasta)

    conn_in, ci = open_input_database(args.database)

    output_chromosomes(genome, ci, basename)