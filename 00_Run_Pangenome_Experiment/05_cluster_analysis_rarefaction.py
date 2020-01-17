#!/usr/bin/env python

'''PGE Rarefaction Summary Plot

Takes the *_binary.tsv file from the 05 PGE Analysis program as input,
calculates the core vs pangenome rarefaction curve with permutations,
and returns a plot as a *_Pangenome_Rarefaction.png file.

This tool takes the following input parameters:

    * bf - binary gene presence absence matrix .tsv file
    * prm - number of permutations to test
    * org - organism name for the plot title (str)
    * op - An output file prefix (str)

This script returns the following files:

    * {op}_plot.pdf
    * {op}_summary.tsv - mean end result of rarefaction after x prms
    * {op}_results.tsv - rarefaction results of each prm

This script requires the following packages:

    * argparse
    * os
    * numpy
    * pandas
    * matplotlib

This file can also be imported as a module and contains the follwing 
functions:

    * pangenome_rarefaction_plot - plots the data. writes the png.
    * gather_data - reads files from input directory. Builds DataFrame
    * main - the main function of the script

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: Thursday, June 12th, 2019
License :: GNU GPLv3
Copyright 2019 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
import os
from collections import defaultdict
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


def Pangenomes_Rarefaction_Plot(df, org, prm, gPg, genomes_per_gene, op):
    ''' Builds Rarefaction Plot and summary then writes them to plot.pdf and summary.txt '''

    oName = ' '.join(org.split('_'))
    x = df.index

    fig, (ax1, ax2) = plt.subplots(
        2, 1,
        gridspec_kw={'height_ratios': [2, 1]},
        figsize=(20, 14),
        sharex=True,
        sharey=False
        )

    ax1.text(
        0.5, 1.13, oName,
        fontstyle='italic', fontweight='heavy', fontsize=60,
        color='#933b41', horizontalalignment='center',
        transform=ax1.transAxes
        )
    ax1.set_title(
        'Core vs Pangenome Rarefaction Curves',
        color='#933b41', fontsize=50, y=1.02
        )
    ax2.set_title(
        'Number of New Genes per Genome / Number of Genes in the Genome',
        color='#a6bddb', fontsize=32, y=1.02
        )

    # Plot Medians, Means and IQR
    ax1.plot(x, df.pan_median, color='#933b41', linestyle='--', lw=1)
    ax1.plot(x, df.pan_mean, color='#933b41', linestyle=':', lw=2)
    ax1.fill_between(x, df.pan025, df.pan975, color='#ddb7b1')

    # Plot Medians, Means and IQR
    ax1.plot(x, df.specific_median, color='#b35806', linestyle='--', lw=1)
    ax1.plot(x, df.specific_mean, color='#b35806', linestyle=':', lw=2)
    ax1.fill_between(x, df.specific025, df.specific975, color='#fdb863')

    # Plot Medians, Means and IQR
    ax1.plot(x, df.rare_median, color='#542788', linestyle='--', lw=1)
    ax1.plot(x, df.rare_mean, color='#542788', linestyle=':', lw=2)
    ax1.fill_between(x, df.rare025, df.rare975, color='#b2abd2')

    # Plot Medians, Means and IQR
    ax1.plot(x, df.common_median, color='#006d2c', linestyle='--', lw=1)
    ax1.plot(x, df.common_mean, color='#006d2c', linestyle=':', lw=2)
    ax1.fill_between(x, df.common025, df.common975, color='#41ab5d')

    # Plot Medians, Means and IQR
    ax1.plot(x, df.core_median, color='#0868ac', linestyle='--', lw=1)
    ax1.plot(x, df.core_mean, color='#0868ac', linestyle=':', lw=2)
    ax1.fill_between(x, df.core025, df.core975, color='#4eb3d3')

    # Set Rarefaction Plot Labels etc.
    t_pan = round(int(list(df.pan_mean)[-1]), 0)
    t_core = round(int(list(df.core_mean)[-1]), 0)
    t_common = round(int(list(df.common_mean)[-1]), 0)
    t_rare = round(int(list(df.rare_mean)[-1]), 0)
    t_specific = round(int(list(df.specific_mean)[-1]), 0)

    stext = (
        f"Pangenome Size: {t_pan}  |  "
        f"Core Genes: {t_core}  |  "
        f"Common Genes: {t_common}  |  "
        f"Rare Genes: {t_rare}  |  "
        f"Specific Genes: {t_specific}"
        )
    ax1.text(
        0.5, -0.05, stext,
        fontsize=18, color='#737373',
        horizontalalignment='center', transform=ax1.transAxes
        )
    ttext = (
        f"Number of Genomes: {len(x)}  |  "
        f"Number of Permutations: {prm}"
        )
    ax1.text(
        0.5, 0.94, ttext,
        fontsize=22, color='#737373',
        horizontalalignment='center', transform=ax1.transAxes
        )

    ax1.set_ylabel('Number of Gene Clusters', fontsize=28)
    ax1.yaxis.grid(which="both", color='#d9d9d9', linestyle='--', linewidth=1)
    ax1.minorticks_on()
    ax1.set_xlabel('')
    ax1.tick_params(labelsize=22)
    for spine in ax1.spines.values(): spine.set_linewidth(2)
    ax1.set_axisbelow(True)

    # Build Rarefaction Plot Legend
    l_pan = Line2D(
        [0],[0], color='w', label='Pan Genome',
        markerfacecolor='#ddb7b1', marker='o', markersize=18
        )
    l_core = Line2D(
        [0],[0], color='w', label='Core Genes',
        markerfacecolor='#4eb3d3', marker='o', markersize=18
        )
    l_common = Line2D(
        [0],[0], color='w', label='Common Genes',
        markerfacecolor='#41ab5d', marker='o', markersize=18
        )
    l_rare = Line2D(
        [0],[0], color='w', label='Rare Genes',
        markerfacecolor='#b2abd2', marker='o', markersize=18
        )
    l_specific = Line2D(
        [0],[0], color='w', label='Specific Genes',
        markerfacecolor='#fdb863', marker='o', markersize=18
        )
    l_mean = Line2D(
        [0],[0], color='#bdbdbd', linestyle=':', lw=4, label='Mean'
        )
    l_median = Line2D(
        [0],[0], color='#bdbdbd', linestyle='--', lw=4, label='Median'
        )
    l_IQR = Line2D(
        [0],[0], color='w', label='95% E.C.I.',
        markerfacecolor='#bdbdbd', marker='s', markersize=20
        )

    legend_elements = [
                        l_pan,
                        l_specific,
                        l_rare,
                        l_common,
                        l_core,
                        l_IQR,
                        l_mean,
                        l_median
                        ]

    ax1.legend(
        handles=legend_elements,
        loc='upper left',
        fontsize=18,
        fancybox=True,
        framealpha=0.0,
        frameon=False
        )

    # Plot Gene Ratio Plot
    ax2.plot(x, df.new_ratio_median, color='#a6bddb', linestyle='--', lw=2)
    ax2.plot(x, df.new_ratio_mean, color='#a6bddb', lw=2)
    ax2.fill_between(x, df.new_ratio025, df.new_ratio975, color='#c9d7e9')

    # Calculate mean new genes per genome without first value
    nRatios = df.new_genes_mean.tolist()[1:]
    mean_new_genes = sum(nRatios) / len(nRatios)
    # Calculate mean new gene ratio without first value of 100%.
    mRatios = df.new_ratio_mean.tolist()[1:]
    mean_new_gene_ratio = sum(mRatios) / len(mRatios)

    mtext = (
        f"Mean Genes per Genome: {int(df.genome_length_mean.mean())}  |  "
        f"Mean New Genes per Genome: {int(mean_new_genes)}"
        )
    ax2.text(
        0.5, 0.90, mtext,
        fontsize=24, color='#737373',
        horizontalalignment='center', transform=ax2.transAxes
        )
    rtext = f'Mean Ratio: {round(mean_new_gene_ratio,2)}%'
    ax2.text(
        len(x)-2, mean_new_gene_ratio+5, rtext,
        fontsize=18, color='#b30000', horizontalalignment='right'
        )
    ax2.plot(
        x, [mean_new_gene_ratio]*len(x),
        linestyle='-.', color='#b30000', lw=1.5
        )
    ax2.set_ylabel('New Genes Ratio', fontsize=28)
    ax2.set_xlabel('Number of Genomes', fontsize=28)
    ax2.yaxis.grid(which="both", color='#d9d9d9', linestyle='--', linewidth=1)
    ax2.minorticks_on()
    ax2.tick_params(labelsize=22)
    ax2.set_xticks(range(0, len(x)+1, 10))
    ax2.set_xlim(-1, len(x)+1)
    for spine in ax2.spines.values(): spine.set_linewidth(2)
    ax2.set_axisbelow(True)

    plt.subplots_adjust(
        left = 0.09,
        right = 0.98,
        bottom = 0.07,
        top = 0.87,
        hspace = 0.2
        )
    plt.savefig(f'{op}_plot.pdf')
    plt.close()

    XX = org.split('_')
    organism = '_'.join(XX[0:2])
    exp = XX[2]

    # Write out summary file
    summary = [
                organism,
                exp,
                t_pan,
                t_core,
                t_common,
                t_rare,
                t_specific,
                df.genome_length_mean.mean(),
                mean_new_genes,
                ]

    s_out = '\t'.join(str(x) for x in summary)
    with open(f'{op}_summary.tsv', 'w') as o: o.write(f'{s_out}\n')

    # Function end


def Pangenome_Rarefaction_Results(df, org, op):
    '''
    Takes rd dictionary and the output prefix org from
    Pangenome_Rarefaction function. Builds DataFrames from Dictionaries
    and calculate mean, median, stdev, iqr. Writes results to
    *_Pangenome_Rarefaction_Results.tsv and returns DataFrame.
    '''
    big_D = pd.DataFrame()
    for k,v in df.items():
        d = pd.DataFrame(v)
        big_D[f'{k}_median'] = d.median().values
        big_D[f'{k}_mean'] = d.mean().values
        big_D[f'{k}025'] = d.quantile(q=0.025).values
        big_D[f'{k}975'] = d.quantile(q=0.975).values

    df = pd.DataFrame(df)

    df.to_csv(f'{op}_results.tsv', sep='\t') 
    return big_D


def Pangenome_Rarefaction_Genes_per_Genome(df, n):
    '''
    Takes the df2 binary and number of genomes n from Pangenome
    Rarefaction function Calculates the number of genes per genome
    and returns a dictionary of genome id as key and number of genes
    as values.
    '''
    genomes_per_gene = df.sum(axis=1).tolist()
    genes_per_genome = df.sum(axis=0)
    genome_id = genes_per_genome.index
    gene_counts = genes_per_genome.values
    gl = {genome_id[i]: gene_counts[i] for i in range(n)}
    return gl, genomes_per_gene


def Pangenome_Rarefaction(df2, org, prm, op):
    '''
    Takes the df2 binary, a dictionary of genome id's and an output file
    prefix permutations Calculations of the core vs pangenome rarefaction
    curves and percent new genes per genome. Returns a Rarefaction
    Results tsv file and the rarefaction plot as a .png file.
    '''

    # Turn off a false positive warning message.
    pd.options.mode.chained_assignment = None
    # https://stackoverflow.com/questions/20625582/how-to-deal-with-settingwithcopywarning-in-pandas

    df2 = pd.read_csv(df2, sep='\t', index_col=0)
    df2 = df2.set_index('Cluster')
    n = len(df2.columns) # Set the number of genomes.
    gl, genomes_per_gene = Pangenome_Rarefaction_Genes_per_Genome(df2, n)

    # initialize a dictionary to keep track of rarefaction bootstrap counts
    rd = {
        'pan': defaultdict(list),
        'core': defaultdict(list),
        'common': defaultdict(list), 
        'rare': defaultdict(list),
        'specific': defaultdict(list),
        'new_ratio': defaultdict(list),
        'new_genes': defaultdict(list),
        'genome_length': defaultdict(list)
        }

    print('Running', prm, 'permutations:')

    # perform this many permutations on the rarefaction curve calculations
    for j in range(prm):
        print('Running Permutation:', j+1)

        # randomly shuffle order of dataframe columns
        #np.random.seed(j*42)
        #genomes = np.random.choice(df2.columns.tolist(), n, replace=False)
        # shuffling the entire dataframe seems to be faster than np.
        df = df2.sample(frac=1, axis=1, replace=False, random_state=j*42)
        #df = df2.sample(n=1000, axis=1, replace=True, random_state=j*42)
        # get the genomes ids in the new order
        genomes = df.columns.tolist()[:n]

        # iterate through the genomes and calculate pangenome stats for each step after the 1st.
        for i in range(0,n):
            c, r, s = 0.90, 0.20, 1/(i+1)
            x = genomes[:i+1] # select genome names for step i
            y = x[-1] # select last genome name for step i
            d = df2[x] # create dataframe the ith genomes
            d['S'] = d.sum(axis=1).values / (i+1) # add new column with sum of rows / x genomes
            pan = len(d[d.S != 0]) # How large is the pangenome?
            core = len(d[d.S >= c]) # How large is the core genome?
            specific = len(d[d.S == s]) # How large is the specific genome?
            common = len(d[(d.S > r) & (d.S < c) & (d.S != s)]) # How large is the common genome?
            rare = len(d[(d.S > s) & (d.S <= r)]) # How large is the rare genome?

            # Calculate new genes per genome / the number of genes in genome
            if i > 0:
                new_genes = pan - rd['pan'][i-1][j]
                new_ratio = (pan - rd['pan'][i-1][j]) / gl[y] * 100
            else:
                new_genes = pan 
                new_ratio = 100.00

            # Track Genome Lengths
            genome_length = gl[y]

            # Update rd dictionary with new calculations.
            rd['pan'][i].append(pan)
            rd['core'][i].append(core)
            rd['common'][i].append(common)
            rd['rare'][i].append(rare)
            rd['specific'][i].append(specific)
            rd['new_ratio'][i].append(new_ratio)
            rd['new_genes'][i].append(new_genes)
            rd['genome_length'][i].append(genome_length)

            # Testing the values are being calculated properly
            #print(d)
            #print(i+1,round(1/(i+1),2),core,common,rare,specific,pan,core+common+rare+specific==pan)

    ''' Testing the nested dictionaries get built properly.
    print(rd)
    for i in range(n):
        for j in range(permutations):
            print(
                rd['pan'][i][j],
                rd['core'][i][j],
                rd['common'][i][j],
                rd['rare'][i][j],
                rd['specific'][i][j],
                rd['new'][i][j]
                )
    '''

    df = Pangenome_Rarefaction_Results(rd, org, op)
    #print(df)
    gPg = list(gl.values()) # pass list of genes per genome 
    Pangenomes_Rarefaction_Plot(df, org, prm, gPg, genomes_per_gene, op)
    print('Script completed successfully.')


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-bf', '--binary_file',
        help='Please specify the binary.tsv input file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-prm', '--number_permutations',
        help='Please specify how many permutations to run!',
        metavar='',
        type=int,
        required=True
        )
    parser.add_argument(
        '-org', '--organism_name',
        help='Organism_name_experiment ex: Escherichia_coli_0042',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-op', '--output_prefix',
        help='What do you want to name the output file?',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # Run this scripts main function
    print('Calculating the rarefaction curve...')

    Pangenome_Rarefaction(
        args['binary_file'],
        args['organism_name'],
        args['number_permutations'],
        args['output_prefix']
        )

if __name__ == "__main__":
    main()
