#!/usr/bin/env python

'''Plot Whole Genome Stats Across Metagenome Samples.

This script reads a directory of *_genome.tsv files from the script:
07e_MagicBlast_CoveragePlus.py
The directory should only contain *_genome.tsv files.
Plots TAD, ANIr, and Relative Abundance for each genome in each sample
along with the mean and median for each sample.

This tool takes the following input parameters:

    * input dir - containing *_genome.tsv files
    * out - output file name

This script returns the following files:

    * Line plots for TAD, ANIr, and Relative Abundance as a .png file

This script requires the following packages:

    * argparse
    * os.listdir
    * os.path.isfile
    * os.path.join
    * collections.defaultdict
    * matplotlib
    * numpy

This file can also be imported as a module and contains the follwing 
functions:

    * gather_data - reads in the files, parses data to dict
    * plot_data - Builds the plots and writes to out_file.png
    * main - the main function of the script

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: Thursday, September 6th, 2019
License :: GNU GPLv3
Copyright 2019 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
from os import listdir
from os.path import isfile, join
from collections import defaultdict
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

def smpl_convert(smpl, x):
    ''' Adjust sample names for the plot '''
#### The smpl needs to be taylored to the file name of the project #####
#### smpl should be the unique metagenome file ID. The time points can
#### have whatever name you would like.
    if x == 0:
        d = {
            'M13': 'Time-Zero',
            'M15': 'One-Week',
            'M16': 'One-Month'
            }
    elif x == 1:
        d = {
            'Time-Zero': 'M13',
            'One-Week': 'M15',
            'One-Month': 'M16'
            }
#### The smpl needs to be taylored to the file name of the project #####

    timepoint = d[smpl]
    return timepoint


def gather_data(gtd):
    ''' Parses data into a dictionaries. Returns data dicts '''

    tads = defaultdict(lambda: defaultdict(list))
    anirs = defaultdict(lambda: defaultdict(list))
    rels = defaultdict(lambda: defaultdict(list))

    tadmn = defaultdict(list)
    animn = defaultdict(list)
    relmn = defaultdict(list)

    file_list = [f for f in listdir(gtd) if isfile(join(gtd, f))]

    for file in file_list:
        with open(f'{gtd}{file}', 'r') as f:
            header = f.readline()
            X = f.readline().rstrip().split('\t')

#### The smpl needs to be taylored to the file name of the project #####
#### sample needs to match keys in smpl_convert function above.
            file_basename = X[0].split('/')[-1].split('_')
            sample = file_basename[1]
            smpl = smpl_convert(sample, 0)
            genome = file_basename[7]
#### The smpl needs to be taylored to the file name of the project #####

            tad = float(X[1])
            anir = float(X[2][:-1])
            relabnd = float(X[3][:-1])

            tads[genome][smpl].append(tad)
            anirs[genome][smpl].append(anir)
            rels[genome][smpl].append(relabnd)

            tadmn[smpl].append(tad)
            animn[smpl].append(anir)
            relmn[smpl].append(relabnd)

    return tads, anirs, rels, tadmn, animn, relmn


def plot_data(tads, anirs, rels, tadmn, animn, relmn, out):
    ''' Plots the data in dict and write plot to out.png '''

#### The smpl needs to be taylored to the file name of the project #####
    plot_order = ['Time-Zero', 'One-Week', 'One-Month']
#### The smpl needs to be taylored to the file name of the project #####

    # Set the colors
    ctad = '#e5f5e0'
    ctadmn = '#006d2c'
    cani = '#efedf5'
    canimn = '#54278f'
    crel = '#fdd0a2'
    crelmn = '#d94801'
    gridM = '#bdbdbd'
    gridm = '#d9d9d9'

    # Build the plot
    fig, (ax1, ax2, ax3) = plt.subplots(
        1, 3,
        figsize=(25,6),
        sharex=True,
        sharey=False
        )

    # Plot title, labels, and text
    ax1.set_title(
        '80% Truncated Average Depth',
        color=ctadmn, fontsize=24, y=1.02
        )
    ax1.set_ylabel(
        '# bps Recruited / Genome Length', color=ctadmn, fontsize=20
        )
    ax1.set_xlabel(
        #'Metagenome Sample Timepoint', color=ctadmn, fontsize=20, y=-0.5
        ''
        )
    ax2.set_title(
        '% Relative Abundance',
        color=crelmn, fontsize=24, y=1.02
        )
    ax2.set_ylabel(
        '# bps Recruited / Total bps * 100', color=crelmn, fontsize=20
        )
    ax2.set_xlabel(
        #'Metagenome Sample Timepoint', color=crelmn, fontsize=20, y=-0.5
        ''
        )
    ax3.set_title(
        'ANIr > 95%',
        color=canimn, fontsize=24, y=1.02
        )
    ax3.set_ylabel(
        'Average Nucleotide Identity of Reads', color=canimn, fontsize=20
        )
    ax3.set_xlabel(
        #'Metagenome Sample Timepoint', color=canimn, fontsize=20, y=-0.5
        ''
        )

    for ax in fig.axes:
        ax.minorticks_on()
        ax.tick_params(
            which='minor', axis='both', left=False, bottom=False
            )
        ax.tick_params(
                    which='major', axis='both',
                    left=False, bottom=True,
                    size=8, width=5, tickdir='in',
                    labelsize=16, zorder=10
                    )
        ax.yaxis.grid(
            which="minor", color=gridm, linestyle='--',
            linewidth=1, alpha=0.6, zorder=1
            )
        ax.yaxis.grid(
            which="major", color=gridM, linestyle='--',
            linewidth=1.5, alpha=0.4, zorder=1
            )
        for spine in ax.spines.values(): spine.set_linewidth(2)
        ax.set_axisbelow(True)
        ax.set_xticklabels(plot_order, rotation=30, ha='right')

        for genome in tads.keys():
            ys = []
            for smpl in plot_order:
                ys.append(tads[genome][smpl])
            ax1.plot(plot_order, ys, color=ctad)
        ax1.plot(
            plot_order,
            [np.mean(tadmn[smpl]) for smpl in plot_order],
            color=ctadmn, linestyle='--', lw=2
            )

        for genome in anirs.keys():
            ys = []
            for smpl in plot_order:
                ys.append(anirs[genome][smpl])
            ax3.plot(plot_order, ys, color=cani)
        ax3.plot(
            plot_order,
            [np.mean(animn[smpl]) for smpl in plot_order],
            color=canimn, linestyle='--', lw=2
            )

        for genome in rels.keys():
            ys = []
            for smpl in plot_order:
                ys.append(rels[genome][smpl])
            ax2.plot(plot_order, ys, color=crel)
        ax2.plot(
            plot_order,
            [np.mean(relmn[smpl]) for smpl in plot_order],
            color=crelmn, linestyle='--', lw=2
            )

    # Build Plot Legend
    iso_single = Line2D(
        [0],[0], color='#525252',
        linestyle='-', lw=4, label='Single Isolate Genome'
        )
    iso_avg = Line2D(
        [0],[0], color='#525252',
        linestyle='--', lw=4, label='Average of Isolates'
        )
    legend_elements = [iso_single, iso_avg]

    ax1.legend(
        handles=legend_elements, loc='lower center',
        fontsize=18, frameon=False
        )
    ax2.legend(
        handles=legend_elements, loc='lower center',
        fontsize=18, frameon=False
        )
    ax3.legend(
        handles=legend_elements, loc='lower center',
        fontsize=18, frameon=False
        )

    # Set plot axis ranges
    alltad = [i for l in tadmn.values() for i in l]
    allrel = [i for l in relmn.values() for i in l]
    ax1.set_ylim(bottom=0, top=max(alltad)+min(alltad))
    ax2.set_ylim(bottom=0, top=max(allrel)+min(allrel))
    ax3.set_ylim(bottom=95, top=100)

    # adjust layout, save, and close
    fig.set_tight_layout(True)
    plt.savefig(f'{out}.png')
    plt.close()

    print(f'Plots written to {out}.png')

    ### Write means to file
    with open(f'{out}.tsv', 'w') as o:
        header = 'Sample\tTAD\tRelA\tANI\n'
        o.write(header)

        mnTAD = [np.mean(tadmn[smpl]) for smpl in plot_order]
        mnRelA = [np.mean(relmn[smpl]) for smpl in plot_order]
        mnANI = [np.mean(animn[smpl]) for smpl in plot_order]

        for i, s in enumerate(plot_order):

            listout = [
                        smpl_convert(s, 1),
                        f'{mnTAD[i]:.4f}',
                        f'{mnRelA[i]:.4f}',
                        f'{mnANI[i]:.4f}'
                        ]

            lineout = '\t'.join(listout)
            
            o.write(f'{lineout}\n')

    print(f'Mean values for each sample written to {out}.tsv')
    print('Seems like the script finished successfully!\n\n')

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-gtd', '--genome_tsv_dir',
        help='Please specify the directory with *_genome.tsv files!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-out', '--output_file',
        help='Please specify the name to us for the output .png file',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # Run this scripts main function
    print('\n\nGenerating List of Genes and Gene Categories by Cluster...')

    (
    tads, anirs, rels, tadmn, animn, relmn
        ) = gather_data(args['genome_tsv_dir'])

    plot_data(
        tads, anirs, rels, tadmn, animn, relmn, args['output_file']
        )


if __name__ == "__main__":
    main()
