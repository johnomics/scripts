#!/usr/bin/python3 -u

import os
import sys
import argparse
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool

import Oceanic.GenomeData as gd
import Oceanic.Chromosome as chrom
import Oceanic.Stats as stats

from Bio import SeqIO

def get_genome_stats(chromosomes):
    genome_stats = stats.Stats("Genome")

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
    chromosomes = {chromosome_name: chrom.Chromosome(chromosome_name, genome) for chromosome_name, in genome.db.fetchall() if chromosome_name != 0}

    for name in chromosomes:
        mapped_blocks += chromosomes[name].mapped_blocks
        mapped_blocks_length += chromosomes[name].mapped_blocks_length
        placed_blocks += chromosomes[name].placed_blocks
        placed_blocks_length += chromosomes[name].placed_blocks_length
        
    print("Map blocks {} length {} of which {} blocks placed, length {}".format(mapped_blocks, mapped_blocks_length, placed_blocks, placed_blocks_length))

    return chromosomes


def write_block(scaffold, start, genome):
    end = genome.blocks[scaffold][start].end
    if start < end:
        seq = genome.sequences[scaffold][start-1:end]
        seq.id = "{}_{}_{}".format(scaffold, start, end)
    else:
        seq = genome.sequences[scaffold][end-1:start].reverse_complement()
        seq.id = "{}_{}_{}".format(scaffold, end, start)

    seq.description = seq.id
    SeqIO.write(seq, genome.revised_fasta, "fasta")

    if genome.revised_db:
        genome.revised_db.execute("insert into scaffold_map values (?,?,?,?,?,?)",
            [0, -1, scaffold, min(start, end), max(start, end), end-start+1])


def reassemble(chromosomes, genome, args):

    print("Reassembly...")

    pool = ThreadPool(args.threads)
    pool.map(lambda x: chromosomes[x].assemble(args), chromosomes.keys())

    get_genome_stats(chromosomes)
    
    assembly = []
    for chromosome in chromosomes:
        assembly += chromosomes[chromosome].write()

    assembled_blocks = 0
    assembled_length = 0
    for blocklist in assembly:
        for block in blocklist:
            assembled_blocks += 1
            assembled_length += block.length
            del genome.blocks[block.scaffold][block.start]
            if not genome.blocks[block.scaffold]:
                del genome.blocks[block.scaffold]

    print("Assembled {} blocks, length {}".format(assembled_blocks, assembled_length))

    deleted_blocks = 0
    deleted_length = 0
    for blocklist in genome.refuse:
        for block in blocklist:
            deleted_blocks += 1
            deleted_length += block.length
            del genome.blocks[block.scaffold][block.start]
            if not genome.blocks[block.scaffold]:
                del genome.blocks[block.scaffold]

    print("Deleted {} blocks, length {}".format(deleted_blocks, deleted_length))

    left_blocks = 0
    left_length = 0
    offcut_blocks = 0
    offcut_length = 0
    for scaffold in genome.blocks:
#        print(scaffold)
#        for start in sorted(genome.blocks[scaffold]):
#            print(genome.blocks[scaffold][start])
#        if scaffold in genome.sequences:
#            for feature in genome.sequences[scaffold].features:
#                if feature.type == 'gene':
#                    print('\t' + repr(feature))
        if scaffold in genome.offcuts:
#            print(scaffold, 'Offcut')
            for start in genome.blocks[scaffold]:
                offcut_blocks += 1
                offcut_length += genome.blocks[scaffold][start].length
        else:
            for start in genome.blocks[scaffold]:
                left_blocks += 1
                left_length += genome.blocks[scaffold][start].length
            
                if genome.revised_fasta:
                    write_block(scaffold, start, genome)
                
                assembly.append([genome.blocks[scaffold][start]])

    if genome.revised_db:
        genome.revised_conn.commit()
    
    if genome.revised_tsv:
        for origscaffold in sorted(genome.origparts):
            for part in genome.origparts[origscaffold]:
                genome.revised_tsv.write('{}\n'.format(repr(part)))

    print("Removed {} unmapped offcuts, length {}".format(offcut_blocks, offcut_length))
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
    print("Blocks:{}\t".format(blocks), stats.genome(scaffold_lengths))


def get_args():
    parser = argparse.ArgumentParser(description='''Output new FASTA file based on linkage map.
    
        -d database
        -f FASTA file
        -g GFF file
        -e errors file
        -t threads
        -r revised
        -a haplomerger
        -o original TSV
        -p prefix''')

    parser.add_argument('-d', '--database', type=str, required=True)
    parser.add_argument('-f', '--fasta', type=str, required=True)
    parser.add_argument('-g', '--gff', type=str, required=True)
    parser.add_argument('-e', '--errors', type=str)
    parser.add_argument('-t', '--threads', type=int, default=1)
    parser.add_argument('-r', '--revised', type=str)
    parser.add_argument('-a', '--haplomerger', type=str)
    parser.add_argument('-o', '--original', type=str)
    parser.add_argument('-p', '--prefix', type=str)

    return parser.parse_args()


if __name__ == '__main__':
    
    args = get_args()

    genome = gd.GenomeData(args)

    chromosomes = load_map(genome)

    assembly = reassemble(chromosomes, genome, args)

    output(assembly)