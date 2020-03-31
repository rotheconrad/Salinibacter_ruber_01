#!/usr/bin/env python

'''PGE Summary - ANI vs Shared Fragment Correlation.

This script reads a directory of PGE Correlations summary tsv files
and builds a scatter plot of ANI values on the x-axis and the ratio
of shared fragments over total fragments on the y-axis. It plots in
gray the results of the total PGE results and in black the result of
the single specified representative experiment.

This tool takes the following input parameters:

    * csd - correlation summary directory (str)
    * org - organism name for the plot title (str)
    * rep - representative experiment number (int, ex: 0042)
    * op - An output file prefix (str)

This script returns the following files:

    * .png file of plot image with {op}_correlation_plot.png

This script requires the following packages:

    * argparse
    * os
    * numpy
    * pandas
    * matplotlib

This file can also be imported as a module and contains the following 
functions:

    * ANI_vs_shared_plot - plots the data. writes the png.
    * gather_data - reads files from input directory. Builds DataFrame
    * gather_stats - computes summary stats and correlation on the data
    * main - the main function of the script

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: Thursday, June 27th, 2019
License :: GNU GPLv3
Copyright 2019 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
import os
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np, pandas as pd; np.random.seed(0)
import seaborn as sns; sns.set(style="white", color_codes=True)


def gather_stats(df):
    """Computes correlation, mean, and median on df columns xs and ys

    Parameters
    ----------
    df : DataFrame
        DataFrame with columns xs: ANI values, ys: ratio shared fragments

    Returns
    -------
    df_stats : dict
        dictionary with keys scorr, pcorr, kcorr, ani_mean, ani_median,
        frag_mean, frag_median
    """

    # Compute Pearson Correlation Coefficient
    pcorr = round(df['xs'].corr(df['ys'], method='pearson'), 3)
    # Compute Spearman's Correlation
    scorr = round(df['xs'].corr(df['ys'], method='spearman'), 3)
    # Compute Kendall's Ta3
    kcorr = round(df['xs'].corr(df['ys'], method='kendall'), 3)
    # Compute ANI mean and median
    ani_mean = np.mean(df['xs'])
    ani_median = np.median(df['xs'])
    # Compute shared fragment mean and median
    frag_mean = np.mean(df['ys'])
    frag_median = np.median(df['ys'])
    # Compile dictionairy
    df_stats = {
        'pcorr': pcorr,
        'scorr': scorr,
        'kcorr': kcorr,
        'ani_mean': ani_mean,
        'ani_median': ani_median,
        'frag_mean': frag_mean,
        'frag_median': frag_median
        }

    print(f"ANI mean: {ani_mean}", f"ANI median: {ani_median}")
    print(f"Frag mean: {frag_mean}", f"Frag median: {frag_median}")

    return df_stats


def ANI_correlation_plot(df_full, org, op):
    """Takes the data and builds the plot

    Parameters
    ----------
    df_full : DataFrame
        Dataframe with data for all the experients run
    df_rep : DataFrame
        Dataframe with data for the selected representative experiment.
    org : str
        The name of the organism for the plot title
    op : str
        The prefix for the output plot.png

    Returns
    -------
    No return. Saves the plot to {op}_correlation_plot.png in the
    working directory unless another directory was specified as
    part of the op.
    """

    org_name = ' '.join(org.split('_'))

    # Gather Stats
    df_full_stats = gather_stats(df_full)
    stats_line = (
        "Pearson Correlation Coefficient:\n"
        "Spearman's Correlation:\n"
        "Kendall's Tau:\n"
        )
    full_stats_line = (
        f"{df_full_stats['pcorr']}\n"
        f"{df_full_stats['scorr']}\n"
        f"{df_full_stats['kcorr']}"
        )

    # Set Colors
    full_color = "#bdbdbd"
    full_line_color = "#000000"
    rep_color = "#000000"

    # Build the plot
    #fig, ax = plt.subplots(figsize=(20, 10))

    # plot the data for df_full
    g = sns.jointplot(
        x=df_full['xs'], y=df_full['ys'],
        marker='o', s=30, color=full_color, alpha=0.15
        )
    full_mn = g.ax_joint.axvline(
            x=df_full_stats['ani_mean'], ymin=0, ymax=1,
            color=full_line_color, linewidth=3, linestyle='--',
            label='Mean'
            )
    full_md = g.ax_joint.axvline(
            x=df_full_stats['ani_median'], ymin=0, ymax=1,
            color=full_line_color, linewidth=3, linestyle=':',
            label='Median'
            )
    _ = g.ax_joint.axhline(
        y=df_full_stats['frag_mean'], xmin=0, xmax=1,
        color=full_line_color, linewidth=3, linestyle='--',
        label='Mean'
        )
    _ = g.ax_joint.axhline(
        y=df_full_stats['frag_median'], xmin=0, xmax=1,
        color=full_line_color, linewidth=3, linestyle=':',
        label='Median'
        )

    # plot title, labels, and text
    g.ax_marg_x.text(
        0.5, 1.53, org_name,
        fontstyle='italic', fontweight='heavy', fontsize=60, color=rep_color,
        horizontalalignment='center', transform=g.ax_marg_x.transAxes
        )
    g.ax_marg_x.set_title(
        'ANI Value vs Shared Sequence Fragment Ratio',
        fontsize=50, y=1.02, color=rep_color
        )
    g.ax_joint.text(
        0.26, 0.99, stats_line,
        fontsize=18, color='#252525',
        verticalalignment='top', horizontalalignment='right',
        transform=g.ax_joint.transAxes
        )
    g.ax_joint.text(
        0.265, 0.99, f'{full_stats_line}',
        fontsize=18, fontweight='bold', color='#252525',
        verticalalignment='top', horizontalalignment='left',
        transform=g.ax_joint.transAxes
        )
    g.ax_joint.set_xlabel(
        'Average Nucleotide Identity (fastANI)',
        fontsize=30, fontweight='bold', y=-0.02
        )
    g.ax_joint.set_ylabel(
        'Shared / Total Fragments',
        fontsize=30, fontweight='bold', x=-0.02
        )

    # set the axis parameters / style
    g.ax_joint.minorticks_on()
    g.ax_joint.set_xticks(np.arange(94.9, 100.1, 0.1), minor=True)
    g.ax_joint.set_xticks(np.arange(95.0, 100.1, 0.5))
    g.ax_joint.set_yticks(np.arange(0.46, 1.02, 0.02), minor=True)
    g.ax_joint.set_yticks(np.arange(0.5, 1.1, 0.1))
    g.ax_joint.set_xlim(left=94.9, right=100.1)
    g.ax_joint.set_ylim(bottom=0.46, top=1.02)
    g.ax_joint.tick_params(axis='both', labelsize=18)
    g.ax_joint.tick_params(
        axis='x', which='major', direction='in', color='k',
        width=6, length=12, bottom=True, zorder=3
        )

    # set grid style
    g.ax_joint.yaxis.grid(which="minor", color='#d9d9d9', linestyle='--', linewidth=1)
    g.ax_joint.xaxis.grid(which="minor", color='#f0f0f0', linestyle='-', linewidth=1)
    g.ax_joint.yaxis.grid(which="major", color='#d9d9d9', linestyle='--', linewidth=1.5)
    g.ax_joint.xaxis.grid(which="major", color='#f0f0f0', linestyle='-', linewidth=2)
    g.ax_joint.set_axisbelow(True)

    g.fig.set_figwidth(20)
    g.fig.set_figheight(10)

    # plot the legend
    g.ax_joint.legend(
        handles=[full_mn, full_md],
        loc='lower right',
        fontsize=24,
        frameon=False,
        ncol=2
        )

    # adjust layout, save, and close
    #g.set_tight_layout(True)
    g.savefig(f'{op}_correlation_plot.png')
    plt.close()    


def gather_data(csd):
    """Gets the data and builds a DataFrame

    Parameters
    ----------
    csd : str
        The directory of PGE Correlation Summary tsv files
    rep : str
        The experiment number to use as the representative data

    Returns
    -------
    df_full : DataFrame
        DataFrame with two columns 'xs' and 'ys' containing all data
    df_rep : DataFrame
        DataFrame with tow columns 'xs' and 'ys' containing only data
        from the specified representative experiment
    """

    data_full = {'xs': [], 'ys': []}

    file_list = os.listdir(csd)

    for file in file_list:

        try:
            with open(f'{csd}{file}', 'r') as f:
                # We do not need the first line here
                _ = f.readline()
                # Second line contains ANI values for x-axis
                xs = f.readline().rstrip().split('\t')
                # Third line contains shared fragment ratio for y-axis
                ys = f.readline().rstrip().split('\t')

                data_full['xs'].extend([float(x) for x in xs])
                data_full['ys'].extend([float(y) for y in ys])

        except:
            print(file)

    df_full = pd.DataFrame(data_full)
    df_full = df_full[df_full['xs'] != 100.0]

    print(df_full.describe())

    return df_full


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-csd', '--correlation_summary_directory',
        help='Please specify the PGE correlation summary directory!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-org', '--organism_name',
        help='Organism name for plot title. ex: Escherichia_coli',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-op', '--output_file_prefix',
        help='What do you want to name the ouput file?',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('Running Script...')
    df_full = gather_data(
        args['correlation_summary_directory'],
        )
    ANI_correlation_plot(
        df_full,
        args['organism_name'],
        args['output_file_prefix']
        )

if __name__ == "__main__":
    main()
