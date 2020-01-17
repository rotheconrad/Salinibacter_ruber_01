#!/usr/bin/env python

'''Takes *.ess.faa files for a collection of genomes from MiGA
and rearranges the sequences into fasta files for each essential gene
that contains the gene sequences for each genome.

Transpose fasta file for each genome containing essential genes to
fasta file for each essential gene containing sequence for each genome.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: Friday, August, 9th, 2019
License :: GNU GPLv3
Copyright 2019 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse, os
from collections import defaultdict
from statistics import median


def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


def collect_essentials(ess):

    genes = defaultdict(list)
    lengths = defaultdict(list)

    ess_dir = ess.split('/')[0] + '/'
    f_list = [f for f in os.listdir(ess_dir) if os.path.isfile(f'{ess_dir}{f}')]

    for file in f_list:
        with open(f'{ess_dir}{file}', 'r') as f:
            for name, seq in read_fasta(f):
                X = name.rstrip().split('|')
                newName = X[0]
                key = X[1]
                entry = f'{newName}\n{seq}\n'
                genes[key].append(entry)
                lengths[key].append(len(seq))

    return genes, lengths


def make_outputdir(opd):
    new_dir = opd.split('/')[0] + '/'
    if not os.path.exists(new_dir) and not os.path.isdir(new_dir):
        os.mkdir(new_dir)
    return new_dir

def transpose_essentials(ess, opd):

    # Check if opd exists and create if not
    out_dir = make_outputdir(opd)

    # Collect all gene sequences and lengths in dict with gene name as key
    genes, lengths = collect_essentials(ess)

    # for each gene write a new file with the sequence of each genome
    for key, value in genes.items():
        with open(f'{out_dir}{key}.ess.faa', 'w') as o:
            for entry in value:
                o.write(entry)

    # write stat summary file with number of genomes that have each gene
    # and the min, mean, median, max gene lengths
    with open('Collect_MiGA_Ess_Summary.tsv', 'w') as o:
        header = (
            'Gene\tNumber of Genomes\tGene Lengths:\tMin\tMean\t'
            'Median\tMax\tDifference (Max-Min)\tPercent Difference\n'
            )

        o.write(header)

        for key, value in lengths.items():
            gene = key
            gsnmb = len(value) 
            gsmin = min(value)
            gsmn = sum(value) / gsnmb
            gsmdn = median(value)
            gsmax = max(value)
            diff = gsmax - gsmin
            diffratio = diff / gsmdn * 100

            line_out = (
                        f'{gene}\t{gsnmb}\t\t{gsmin}\t{gsmn:.2f}'
                        f'\t{gsmdn}\t{gsmax}\t{diff}\t{diffratio:.2f}\n'
                        )

            o.write(line_out)

    #fin


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-ess', '--ess.faa_dir',
        help='Directory containing the *ess.faa files to transpose.',
        metavar='',
        type=str,
        #required=True
        )
    parser.add_argument(
        '-opd', '--out_put_directory',
        help='Directory to write new files to.',
        metavar='',
        type=str,
        #required=True
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('Running Script...')
    transpose_essentials(
                        args['ess.faa_dir'],
                        args['out_put_directory']
                        )


if __name__ == "__main__":
    main()
