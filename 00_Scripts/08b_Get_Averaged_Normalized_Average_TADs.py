#!/usr/bin/env python

'''Get Averaged Normalized Average TAD for each gene cluster.

Part 1
Reads the TADgene directory from the 07e_MagicBlast_CoveragePlus.py
script and matches the TAD and gene length values to the
00_Genes_Clusters_PanCats.tsv file from the
08a_Get_Genes_Clusers_PanCat.py script.

Part 2
Calculate Average TAD for each Gene Cluster

Part 3
Normalize averaged gene cluster TAD by average genome TAD of each sample.

Part 4
Average the averaged normalized TAD of each gene cluster across samples.

This tool takes the following input parameters:

    * TADgene/ - directory containing *gene_tad.tsv files from
      07e_MagicBlast_CoveragePlus.py
    * Genes_Clusters_PanCats.tsv from 08a_Get_Genes_Clusters_PanCat.py
    * out - output file name

This script returns the following files:

    * 01_output_file_name_gene_tad_pancat.tsv
    * 02_output_file_name_Part_2-4.tsv
    * 03_output_file_name_Part_4.png

This script requires the following packages:

    * argparse
    * os.makedirs
    * os.listdir
    * os.path.exists
    * os.path.isfile
    * os.path.join
    * collections.defaultdict
    * matplotlib

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: January 21st, 2020
License :: GNU GPLv3
Copyright 2020 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
from os import makedirs, listdir
from os.path import exists, isfile, join
from collections import defaultdict
#import matplotlib
#matplotlib.use('agg')
#import matplotlib.pyplot as plt

def get_gene_tad_pancat(taddir, pancatfile, outdir):
    ''' Gets TAD and gene legnth for each gene of each genome.
        Combines with Gene Cluster and Pancat info.
        Writes tsv file: Gene_Cluster_PanCat_n/N_GeneLen_mtgxTAD.
        Returns dict of {Gene: list} where list is the tsv data '''

    metagenome_tad = defaultdict(lambda: defaultdict())
    gene_lengths = {}
    gene_tad_pancat = {}

    file_list = [f for f in listdir(taddir) if isfile(join(taddir, f))]

    for file in file_list:
        file_base_name = file.split('/')[-1]
        mtg = '_'.join(file_base_name.split('_')[:2])

        with open(f'{taddir}{file}', 'r') as f:
            headerA = f.readline()
            for l in f:
                X = l.rstrip().split('\t')
                gene = X[0]
                tad = X[1]
                glen = X[2]
                metagenome_tad[mtg][gene] = tad
                gene_lengths[gene] = glen

    outfilename = f'{outdir}01_PanCats_ANATADs_byGene.tsv'

    with open(pancatfile, 'r') as f, open(outfilename, 'w') as o:
        header = f.readline().rstrip().split('\t')
        header.append('Gene-Length')
        for m in metagenome_tad.keys():
            header.append(f'{m}-TAD')

        o.write('\t'.join(header) + '\n')

        for l in f:
            X = l.rstrip().split('\t')
            gene = X[0]
            X.append(gene_lengths[gene])

            for metagenome, TADS in metagenome_tad.items():
                X.append(TADS[gene])

            gene_tad_pancat[gene] = X
            o.write('\t'.join(X) + '\n')

    return gene_tad_pancat


def get_cluster_averages(data, outdir):
    ''' Averages the TAD by cluster for each metagenome (sample) '''

    # dictionary for cluster tads for each metagenome.
    # {metagenome: {cluster: [tads]}}
    clstr_tads = defaultdict(lambda: defaultdict(list))
    clstr_lens = defaultdict(list) # dictionary for cluster lens
    clstr_pancat = {} # dictionary for cluster category and value

    for gene, entry in data.items():
        cluster = entry[1]
        pancat = entry[2:4]
        glen = entry[4]
        tads = entry[5:]

        for i, t in enumerate(tads):
            mtg = f'mtg-0{i}'
            print(mtg, t)
            clstr_tads[mtg][cluster].append[t]

        clstr_lens[cluster].append(glen)
        clstr_pancat[cluster] = pancat

        print(cluster, pancat, glen, tads)

def normalize_average_metagenomes(cluster_avgtad, gStats):
    pass


def write_cluster_data(data):
    pass


def plot_cluster_data(data):
    pass


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-t', '--TADgene_dir',
        help='Please specify the directory with *_gene_tad.tsv files!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-g', '--Genes_Clusters_PanCats',
        help='Output file from 08a_Get_Genes_Clusters_PanCat.py!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-s', '--Genome_Stats',
        help='Output tsv from 07h_Genome_Stats.py!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_directory',
        help='Please specify the output directory!',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # Run this scripts main function
    print('\n\nRunning Script ...\n')

    # Check output directory and create it if needed.
    outdir = args['output_directory']
    if not outdir[-1] == '/': outdir = f'{outdir}/'
    if not exists(outdir): makedirs(outdir)

    # PART 1
    print('01: Collecting Gene TADs and matching to PanCats.\n')
    gene_tad_pancat = get_gene_tad_pancat(
                                        args['TADgene_dir'],
                                        args['Genes_Clusters_PanCats'],
                                        outdir
                                        )

    # PART 2
    print('02: Average TADs per cluster for each metagenome (sample).\n')
    cluster_avgtad = get_cluster_averages(gene_tad_pancat, outdir)

    # PART 3
    print('03: Normalizing cluster TADs by metagenome and averaging samples.\n')
    normalized = normalize_average_metagenomes(
                                        cluster_avgtad,
                                        args['Genome_Stats']
                                        )

    # PART 4 - write and plot data from parts 2 & 3.
    print('04a: Writing Normalized, Averaged, TADs by Gene Cluster.')
    _ = write_cluster_data(normalized)

    print('04: Building Summary Plot..\n\n')
    _ = plot_cluster_data(normalized)

if __name__ == "__main__":
    main()
