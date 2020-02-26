#!/usr/bin/env python

''' Plot Gene TAD changes across metagenome timepoints.

This script takes four inputs:
    * Output tsv file from 08c_Match_ANA-TADs_toAnnotations.py script
    * txt file containing column names to plot in order (1 per line).
      These are the names of the columns with the "-NORM" in them in the
      order you would like them plotted. This file should contain two
      columns separated by a comma:
      sampleName-NORM, Label name for x-axis
      with one *-NORM name per line.
    * txt file containing cluster numbers to plot(1 per line). This file
      should contain two columns separated by a comma:
      cluster_number, hexidecimal color(#FFFFFF)
      with one cluster_number per line.
      Checkout for http://colorbrewer2.org/ for color scheme ideas.
    * Category of genes to plot: (Specific, Rare, Common, Core, All)

This script returns a line plot in .png format of the cluster numbers
provided plotted in color and the remaining clusters from the
user-specified pangenome category (Specific, Rare, Common, Core, All)
plotted in light grey.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: February 24th, 2020
License :: GNU GPLv3
Copyright 2020 Roth Conrad
All rights reserved
-------------------------------------------
'''


import argparse
from collections import defaultdict
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt


def read_input_file(infile, clusters, morder, pancat):
    '''Reads input file and returns dictionary with data for plot'''

    # intitialize structures
    main_labels = [0 for i in clusters]
    main_values = [0 for i in clusters]
    other_values = []
    corder = []

    # select PanCat
    if pancat == 'All': pancat = ['Rare', 'Common', 'Core', 'Specific']
    else: pancat = [pancat]

    with open(infile, 'r') as f:
        header = f.readline().rstrip().split('\t')

        # Get index of TrEMBL annotation column from header.
        trembl = header.index('TrEMBL')

        # Get index numbers in order of morder from header.
        for x in morder:
            corder.append(header.index(x))

        # Read through the file
        for l in f:
            X = l.rstrip().split('\t')
            clust = X[1] # select cluster number
            pc = X[2] # select pancat
            name = X[trembl] # select TrEMBL annotation for legend
            ys = [float(X[i]) for i in corder] # Select *-NORM values in morder

            if clust in clusters:
                pos = clusters.index(clust)
                main_labels[pos] = name
                main_values[pos] = ys
            elif pc in pancat:
                other_values.append(ys)

    return main_labels, main_values, other_values


def read_cluster_file(infile):
    '''Reads the cluster file into two separate lists and returns them'''

    clusters = []
    colors = []

    with open(infile, 'r') as f:
        for l in f:
            X = l.rstrip().split(', ')
            cluster = X[0]
            color = X[1]
            clusters.append(cluster)
            colors.append(color)

    return clusters, colors


def read_metagenome_order(infile):
    '''Reads metagenome order into a list and returns it'''

    morder = []
    xlabels = []

    with open(infile, 'r') as f:
        for l in f:
            X = l.rstrip().split(', ')
            sample = X[0]
            label = X[1]
            morder.append(sample)
            xlabels.append(label)

    return morder, xlabels


def plot_gene_changes(
                    main_labels,
                    main_values,
                    other_values,
                    colors,
                    morder,
                    xlabels,
                    pancat,
                    outfile,
                    ptitle,
                    ymax
                    ):
    ''' This function builds a line plot gene change across metagenomes'''

    # Set the colors
    H1 = '#252525'
    grey = '#d9d9d9'

    # Build the plot
    fig, ax = plt.subplots(figsize=(20, 10))

    # plot title, labels, and text
    ax.set_title(
        ' '.join(ptitle),
        color=H1, fontsize=50, y=1.02
        )
    ax.set_ylabel('Gene TAD / Avg Genome TAD', fontsize=28)
    #ax.set_xlabel('Metagenome Sample / Time Point', fontsize=28)

    # Plot List of Genes with colors
    for i, anno in enumerate(main_labels):
        ax.plot(
            xlabels, main_values[i],
            color=colors[i],
            linestyle='-',
            lw=6,
            label=anno,
            zorder=2
            )

    # Plot Remaining gene category in light grey.
    if pancat == 'All': other = 'All Other Genes'
    else: other = f'Other {pancat} Genes'

    for ys in other_values:
        ax.plot(
            xlabels,
            ys,
            color=grey,
            linestyle='-',
            lw=1,
            label=other,
            zorder=1
            )

    # set the axis parameters / style
    ax.minorticks_on()
    ax.set_yticks(range(0, ymax+1))
    ax.tick_params(axis='both', labelsize=18)
    ax.tick_params(axis='x', labelrotation=15)
    # set grid style
    ax.yaxis.grid(
        which="minor", color='#f0f0f0', linestyle='--', linewidth=1.5
        )
    ax.yaxis.grid(
        which="major", color='#d9d9d9', linestyle='--', linewidth=2
        )
    ax.set_axisbelow(True)

    # Setup the Legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(
        handles[:len(main_labels)+1],
        labels[:len(main_labels)+1],
        markerscale=2,
        ncol=1,
        loc='center left',
        bbox_to_anchor=(1.05, 0.5),
        fancybox=True,
        shadow=True,
        fontsize=18
        )

    # adjust layout, save, and close
    fig.set_tight_layout(True)
    plt.savefig(outfile)
    plt.close()



def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_file',
        help='Please specify the input file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-c', '--clusters_to_plot',
        help='Please specify a file with clusters and colors to plot.',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-m', '--metagenome_order',
        help='Please specify a file with order of *-NORM column names.',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-p', '--pangenome_category',
        help=(
            "Please choose a category of genes to plot. "
            "One of: (Specific, Rare, Common, Core, All."
            ),
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_file_name',
        help='What do you want to name the output file?',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-t', '--plot_title',
        help='What do you want the plot title to be?',
        nargs='+',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-y', '--yaxis_max',
        help='Set max value of y-axis (eg: 9)',
        metavar='',
        type=int,
        required=True
        )
    args=vars(parser.parse_args())

    # Do the damn thing!
    print('\n\nGenerating the plot ...\n')

    # Set variables!
    infile = args['input_file']
    pancat = args['pangenome_category']
    outfile = args['output_file_name']
    ptitle = args['plot_title']
    ymax = args['yaxis_max']

    # Load files and parse data for the plot!
    clusters, colors = read_cluster_file(args['clusters_to_plot'])
    morder, xlabels = read_metagenome_order(args['metagenome_order'])
    main_labels, main_values, other_values = read_input_file(
                                                            infile,
                                                            clusters,
                                                            morder,
                                                            pancat,
                                                            )


    # Plot it up son!
    plot_gene_changes(
                    main_labels,
                    main_values,
                    other_values,
                    colors,
                    morder,
                    xlabels,
                    pancat,
                    outfile,
                    ptitle,
                    ymax
                    )

    print(
        '\n\nCongratulations!! '
        'The script seems to have finished successfully.\n\n'
        )


if __name__ == "__main__":
    main()