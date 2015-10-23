#!/usr/bin/env python3

import os
import glob
import argparse
import statistics

lengths = { '01':17206585, '02':9045316, '03':10716509, '04':9793051, '05':9488091, '06':14090625, '07':14377450, '08':9294343, '09':8872268, '10':17966761, '11':11763673, '12':16327298, '13':18445720, '14':9254488, '15':10474703, '16':10083215, '17':14854129, '18':16816904, '19':16414807, '20':14979251, '21':13824031 }
chromosomes = sorted(lengths.keys())

parser=argparse.ArgumentParser(description='''Sex individuals by coverage of sex chromosome in BAM
    -g bamglob''')
parser.add_argument('-g', '--bamglob', type=str, required=True, default="*.bam")
args = parser.parse_args()

bamfiles = glob.glob(args.bamglob)
for bamfile in sorted(bamfiles):
    print(bamfile + '\t', end='')
    chromosome_reads = {}
    for line in os.popen('samtools view ' + bamfile):
        scaffold = line.split('\t')[2]
        if scaffold.startswith('Hmel2'):
            chromosome = scaffold[5:7]
            if chromosome not in lengths:
                continue
            if chromosome not in chromosome_reads:
                chromosome_reads[chromosome] = 0
            chromosome_reads[chromosome] += 1

    autosome_ratios = []
    for c in chromosomes:
        ratio = chromosome_reads[c] / (lengths[c] / 1000)
        print('{}\t{:.2f}\t'.format(chromosome_reads[c], ratio),end='')
        if c != '21':
            autosome_ratios.append(ratio)
        else:
            mean_aratio = statistics.mean(autosome_ratios)
            print('{:.2f}\t'.format(mean_aratio),end='')
            sex_auto = ratio / mean_aratio
            print('{:.2f}\t'.format(sex_auto))
