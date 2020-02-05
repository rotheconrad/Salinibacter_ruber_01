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

    * 01_PanCats_TADs_byGene.tsv
    * 02_PanCats_ANATADs_byCluster.tsv
    * 03_PanCats_ANATADS_Plots.png

This script requires the following packages:

    * argparse
    * from os import makedirs, listdir
    * from os.path import exists, isfile, join
    * from collections import defaultdict
    * pandas as pd
    * matplotlib (pyplot as plt)
    * seaborn as sns

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
import pandas as pd
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="whitegrid")

# PART 1
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
        #mtg = '_'.join(file_base_name.split('_')[:2])
        mtg = file_base_name.split('_')[1]

        with open(f'{taddir}{file}', 'r') as f:
            headerA = f.readline()
            for l in f:
                X = l.rstrip().split('\t')
                gene = X[0]
                tad = X[1]
                glen = X[2]
                metagenome_tad[mtg][gene] = tad
                gene_lengths[gene] = glen

    outfilename = outdir + '01_PanCats_TADs_byGene.tsv'
    header = '' # define semi global variable outside below local
    smpl_order = [] # Keep metagenome samples in order.

    with open(pancatfile, 'r') as f, open(outfilename, 'w') as o:
        header = f.readline().rstrip().split('\t')
        header.append('Gene-Length')
        for m in metagenome_tad.keys():
            header.append(f'{m}-TAD')
            smpl_order.append(m)

        o.write('\t'.join(header) + '\n')

        for l in f:
            X = l.rstrip().split('\t')
            gene = X[0]
            X.append(gene_lengths[gene])

            for smpl in smpl_order:
                X.append(metagenome_tad[smpl][gene])

            gene_tad_pancat[gene] = X
            o.write('\t'.join(X) + '\n')

    return gene_tad_pancat, header, smpl_order


# PART 2
def get_cluster_averages(data):
    ''' Averages the TAD by cluster for each metagenome (sample) '''

    # dictionary for cluster tads for each metagenome.
    # {metagenome: {cluster: [tads]}}
    clstr_tads = defaultdict(lambda: defaultdict(list))
    clstr_lens = defaultdict(list) # dictionary for cluster lens
    clstr_pancat = {} # dictionary for cluster category and value

    # Read the input dictionary and populate new dictionaries with data
    # organized by gene cluster rather than gene name.
    for gene, entry in data.items():
        cluster = entry[1]
        pancat = entry[1:4]
        glen = int(entry[4])
        tads = entry[5:]

        for i, t in enumerate(tads):
            mtg = f'mtg-0{i}'
            clstr_tads[mtg][cluster].append(float(t))

        clstr_lens[cluster].append(glen)
        clstr_pancat[cluster] = pancat

    # Calculate average values for each gene cluster and metagenome
    # populate new dictionary organized by gene cluster containing
    # combined averages for gene length and TAD for each metagenome
    # columns: Cluster, PanCat, n/N, AvgLen, mtgx-AvgTAD
    clstr_data = {}

    # for each cluster, pancat
    for c, p in clstr_pancat.items():
        # retrieve and calculate average gene length of cluster
        avg_len = f'{sum(clstr_lens[c]) / len(clstr_lens[c]):.2f}'
        # append average gene length to pancat info
        p.append(avg_len)
        # for each metagenome
        for smpl, clust in clstr_tads.items():
            # retrieve and calculate average tad of cluster
            avg_tad = f'{sum(clust[c]) / len(clust[c]):.2f}'
            # append average tad of cluster for metagenome
            p.append(avg_tad)

        # add new entry to clstr_data
        clstr_data[c] = p

    return clstr_data


# PART 3
def normalize_average_metagenomes(data, gStats, header, smpl_order, outdir):
    # keep track of header: add mtgx-Normed
    header = header[1:]
    #avg_genome_tad = []
    avg_genome_tad = {}
    normalized = {}

    # Retrieve the average genome tad for each metagenome
    with open(gStats, 'r') as f:
        ignore_header = f.readline()
        for l in f:
            X = l.split('\t')
            mtg = X[0]
            tad = float(X[1])
            #header.append(f'{mtg}-NORM')
            #avg_genome_tad.append(tad)
            avg_genome_tad[mtg] = tad

    # Append new column names to header:
    for mtg in smpl_order: header.append(f'{mtg}-NORM')
    # append column name to header
    # averaged, normalized, average tad (ANA-TAD)
    # averaged gene cluster tad per metagenome
    # normalized by average genome tad per metagenome
    # averaged across metagenomes
    header.append('ANA-TAD')

    # divide tad for each cluster for each metagenome by avg_genome_tad
    for clust, entry in data.items():
        tads = entry[4:]
        norms = []
        for i, mtg in enumerate(smpl_order):
            normed = float(tads[i]) / avg_genome_tad[mtg]
            norms.append(normed)
        # average the averages
        anatad = sum(norms) / len(norms)
        # append new columns to data: mtgx-Norm, ANA-TAD
        for n in norms: entry.append(f'{n:.4f}')
        entry.append(f'{anatad:.4f}')
        # populate dictionary with new values {cluster: data}
        normalized[clust] = entry

    # Write new file for cluster ANATADs
    outfilename = outdir + '02_PanCats_ANATADs_byCluster.tsv'
    with open(outfilename, 'w') as o:
        o.write('\t'.join(header) + '\n')
        for entry in normalized.values():
            o.write('\t'.join(entry) + '\n')

    return outfilename


# PART 4
def plot_cluster_data(df, outdir):
    ''' Plots the data writes png to output_file '''

    # Set the colors
    ccore = '#4eb3d3'
    ccommon = '#41ab5d'
    crare = '#b2abd2'
    cspecific = '#fdb863'
    colors = [ccore, ccommon, crare, cspecific]
    extended_colors = [cspecific] + ([crare]*19) + ([ccommon]*70) + ([ccore]*10)
    gridminor = '#f0f0f0'
    gridmajor = '#d9d9d9'

    # Build the plot
    fig, (
        ax1, ax2, ax3
        ) = plt.subplots(
                        3, 1,
                        figsize=(45, 25),
                        sharex=True,
                        sharey=False
                        )

    # Extended boxplot for each n/N
    ax1 = sns.boxplot(
        x="n/N", y="ANA-TAD", data=df,
        palette=extended_colors, ax=ax1
        )
    ax1.set_xlabel(
        ''
        )
    ax1.set_ylabel(
        'ANA-TAD of Gene Cluster',
        fontsize=24, fontweight='bold', x=-0.02
        )
    ax2 = sns.countplot(
        x="n/N", data=df,
        palette=extended_colors, ax=ax2
        )
    ax2.set_xlabel(
        ''
        )
    ax2.set_ylabel(
        'Number of Gene Clusters in n/N Fraction',
        fontsize=24, fontweight='bold', x=-0.02
        )
    ax3 = sns.boxplot(
        x="n/N", y="Gene-Length", data=df,
        palette=extended_colors, ax=ax3
        )
    xlabel = (
        'Fraction of Genomes with Gene Cluster (n/N) where:\n'
        'n/N = number of genomes with gene in cluster / total genomes'
        )
    ax3.set_xlabel(
        xlabel,
        fontsize=48, fontweight='bold', y=-0.1
        )
    ax3.set_ylabel(
        'Average Gene Length of Gene Cluster',
        fontsize=24, fontweight='bold', x=-0.02
        )

    # set the axis parameters / style
    for ax in fig.axes:
        ax.minorticks_on()
        ax.tick_params(axis='both', labelsize=18)
        ax.tick_params(axis='x', labelrotation=45)
        # set grid style
        ax.yaxis.grid(
            which="minor", color=gridminor, linestyle='--', linewidth=2.5
            )
        ax.yaxis.grid(
            which="major", color=gridmajor, linestyle='--', linewidth=3
            )
        ax.set_axisbelow(True)

    #plt.xticks(rotation=45)
    # adjust layout, save, and close
    fig.set_tight_layout(True)
    outfilename = outdir + '03_PanCats_ANATADs_Plots.png'
    plt.savefig(outfilename)
    plt.close() 


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
    gene_tad_pancat, header, smpl_order = get_gene_tad_pancat(
                                        args['TADgene_dir'],
                                        args['Genes_Clusters_PanCats'],
                                        outdir
                                        )

    # PART 2
    print('02: Average TADs per cluster for each metagenome (sample).\n')
    cluster_avgtad = get_cluster_averages(gene_tad_pancat)

    # PART 3
    print('03: Normalizing cluster TADs by metagenome and averaging samples.\n')
    outfilename = normalize_average_metagenomes(
                                        cluster_avgtad,
                                        args['Genome_Stats'],
                                        header,
                                        smpl_order,
                                        outdir
                                        )

    # PART 4 - Build some plots.
    print('04: Building Summary Plots ...\n\n')
    df = pd.read_csv(outfilename, sep='\t')
    _ = plot_cluster_data(df, outdir)


if __name__ == "__main__":
    main()
