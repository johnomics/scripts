#!/usr/bin/env python3

import os
import glob
import argparse
import statistics

lengths = { '01':17206585, '02':9045316, '03':10716509, '04':9793051, '05':9488091, '06':14090625, '07':14377450, '08':9294343, '09':8872268, '10':17966761, '11':11763673, '12':16327298, '13':18445720, '14':9254488, '15':10474703, '16':10101442, '17':14854129, '18':16816904, '19':16414807, '20':14979251, '21':13824031 }
chromosomes = sorted(lengths.keys())

parser=argparse.ArgumentParser(description='''Sex individuals by coverage of sex chromosome in BAM
    -g bamglob
    -a agp
    -s sex_scaffold_output''')
parser.add_argument('-g', '--bamglob', type=str, required=True, default="*.bam")
parser.add_argument('-a', '--agp', type=str, required=False)
parser.add_argument('-s', '--sex_scaffold_output', action='store_true')
args = parser.parse_args()

sex_scaffold_lengths = {}
if args.agp:
    for line in open(args.agp, 'r'):
        f = line.rstrip().split('\t')
        if 'unmapped' in f[0] or f[4] != 'D':
            continue
        chromosome, start, end, part, parttype, scaffold, scaffold_start, scaffold_end, *fargs = f
        chromosome = chromosome[3:]
        if len(chromosome) == 1:
            chromosome = '0' + chromosome
        lengths[chromosome] = int(end)
        if scaffold.startswith('Hmel221'):
            sex_scaffold_lengths[scaffold] = int(scaffold_end)
    if args.sex_scaffold_output:
        print('File',end='')
        for s in sorted(sex_scaffold_lengths):
            print('\t{}'.format(s),end='')
        print()
        print('Length',end='')
        for s in sorted(sex_scaffold_lengths):
            print('\t{}'.format(sex_scaffold_lengths[s]),end='')
        print()

bamfiles = glob.glob(args.bamglob)
for bamfile in sorted(bamfiles):
    output = bamfile + '\t'
    chromosome_reads = {}
    sex_scaffold_reads = {}
    for line in os.popen('samtools view ' + bamfile):
        scaffold = line.split('\t')[2]
        if scaffold.startswith('Hmel2'):
            chromosome = scaffold[5:7]
            if chromosome not in lengths:
                continue
            if chromosome not in chromosome_reads:
                chromosome_reads[chromosome] = 0
            chromosome_reads[chromosome] += 1
            
            if scaffold.startswith('Hmel221'):
                if scaffold not in sex_scaffold_reads:
                    sex_scaffold_reads[scaffold] = 0
                sex_scaffold_reads[scaffold] += 1


    autosome_ratios = []
    sex_auto = mean_aratio = 0
    for c in chromosomes:
        ratio = chromosome_reads[c] / (lengths[c] / 1000)
        output += '{}\t{:.2f}\t'.format(chromosome_reads[c], ratio)
        if c != '21':
            autosome_ratios.append(ratio)
        else:
            mean_aratio = statistics.mean(autosome_ratios)
            output += '{:.2f}\t'.format(mean_aratio)
            sex_auto = ratio / mean_aratio
            output += '{:.2f}'.format(sex_auto)
    
    if args.sex_scaffold_output:
        output = bamfile
        for s in sorted(sex_scaffold_lengths):
            if s not in sex_scaffold_reads:
                sex_scaffold_reads[s] = 0
            ratio = (sex_scaffold_reads[s] / (sex_scaffold_lengths[s] / 1000)) / mean_aratio
            output += '\t{:.2f}'.format(ratio)
    
    print(output)
