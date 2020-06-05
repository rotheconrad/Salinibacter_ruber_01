#!/usr/bin/env python

'''MiGA All vs All ANI vs Shared Fragment Correlation.

This script reads the all vs all ANI file from MiGA, computes the
mean, meadian and correlation of ANI values (x-axis) and Shared / Total
Fragments (y-axis). Returns a scatter plot of the data as a png file.

This tool takes the following input parameters:

    * ani - all vs all ani file from MiGA (str)
    * org - organism name for the plot title (str)
    * op - An output file prefix (str)

This script returns the following files:

    * .png file of plot image with {op}_correlation.pdf
    * .tsv file with 3 lines.
        #Line 1: Organism Pearson Spearman Kendall Correlation
        #Line 2: xs (ANI values)
        #Line 3: ys (Shared / Total Fragments)

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
Date Created :: Thursday, June 21st, 2019
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


def gather_stats(df, org, op):
    """Computes correlation, mean, and median on df columns xs and ys

    Parameters
    ----------
    df : DataFrame
        DataFrame with columns xs: ANI values, ys: ratio shared fragments
    org : str
        The name of the organism for the plot title
    op : str
        The prefix for the output plot.png

    Returns
    -------
    df_stats : dict
        dictionary with keys scorr, pcorr, kcorr, ani_mean, ani_median,
        frag_mean, frag_median
    {op}_summary.tsv : file
        *.tsv file with 3 lines.
        #Line 1: Organism Pearson Spearman Kendall Correlation
        #Line 2: xs (ANI values)
        #Line 3: ys (Shared / Total Fragments)
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
    '''
    with open(f'{op}_summary.tsv', 'w') as o:
        header=f'{org}\t{pcorr}\t{scorr}\t{kcorr}\n'
        xline='\t'.join([str(x) for x in df['xs']])
        yline='\t'.join([str(y) for y in df['ys']])
        o.write(header)
        o.write(f'{xline}\n')
        o.write(f'{yline}\n')
    '''
    return df_stats


def ANI_vs_shared_plot(df, org, op):
    """Takes the data and builds the plot

    Parameters
    ----------
    df : DataFrame
        Dataframe with data for the experients
    org : str
        The name of the organism for the plot title
    op : str
        The prefix for the output plot.png

    Returns
    -------
    No return. Writes two files.
        #File 1: op + _summary.tsv
        #File 2: op + _plot.png
    """

    org_name = ' '.join(org.split('_'))

    # Gather Stats
    df_stats = gather_stats(df, org, op)
    stats_line1 = (
        "Pearson Correlation Coefficient:\n"
        "Spearman's Correlation:\n"
        "Kendall's Tau:\n"
        )
    stats_line2 = (
        f"{df_stats['pcorr']}\n"
        f"{df_stats['scorr']}\n"
        f"{df_stats['kcorr']}"
        )

    # Set Colors
    main_color = '#933b41'
    second_color = '#737373'
    vline_color = '#000000'

    # Build the plot

    # Plot data points
    '''
    g = sns.jointplot(
        x=df1['xs'], y=df1['ys'],
        marker='o', s=30, color=main_color, alpha=0.15
        )
    '''
    colors = ['#66c2a5','#8da0cb','#fc8d62','#e78ac3','#a6d854','#ffd92f']

    g = sns.JointGrid("xs", "ys", df)

    for i, (cat, Comp) in enumerate(df.groupby("Comp")):

        sns.kdeplot(
                Comp["xs"],
                ax=g.ax_marg_x,
                legend=False,
                color=colors[i]
                )
        sns.kdeplot(
                Comp["ys"],
                ax=g.ax_marg_y,
                vertical=True,
                legend=False,
                color=colors[i]
                )
        g.ax_joint.plot(
                Comp["xs"],
                Comp["ys"],
                "o",
                ms=5,
                alpha=0.15,
                color=colors[i],
                label=cat
                )



    # plot title, labels, text    
    g.ax_marg_x.text(
        0.5, 1.53, org_name,
        fontstyle='italic', fontweight='heavy', fontsize=60, color=main_color,
        horizontalalignment='center', transform=g.ax_marg_x.transAxes
        )
    g.ax_marg_x.set_title(
        'ANI Value vs Shared Sequence Fragment Ratio',
        fontsize=50, y=1.02, color=main_color
        )
    g.ax_joint.text(
        0.26, 0.99, stats_line1,
        fontsize=18, color=second_color,
        verticalalignment='top', horizontalalignment='right',
        transform=g.ax_joint.transAxes
        )
    g.ax_joint.text(
        0.265, 0.99, stats_line2,
        fontsize=18, fontweight='bold', color=main_color,
        verticalalignment='top', horizontalalignment='left',
        transform=g.ax_joint.transAxes
        )
    g.ax_joint.set_xlabel(
        'Average Nucleotide Identity (MiGA ANI)',
        fontsize=30, fontweight='bold', y=-0.02
        )
    g.ax_joint.set_ylabel(
        'Shared / Total Fragments',
        fontsize=30, fontweight='bold', x=-0.02
        )

    # set the axis parameters / style
    g.ax_joint.minorticks_on()
    g.ax_joint.set_xticks(np.arange(93.9, 100.1, 0.1), minor=True)
    g.ax_joint.set_xticks(np.arange(94.0, 100.1, 0.5))
    g.ax_joint.set_yticks(np.arange(0.36, 1.02, 0.02), minor=True)
    g.ax_joint.set_yticks(np.arange(0.5, 1.1, 0.1))
    g.ax_joint.set_xlim(left=93.9, right=100.1)
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


    # Plot mean and median
    mn = g.ax_joint.axvline(
        x=df_stats['ani_mean'], ymin=0, ymax=1,
        color=vline_color, linewidth=2, linestyle='--',
        label='Mean'
        )
    md = g.ax_joint.axvline(
        x=df_stats['ani_median'], ymin=0, ymax=1,
        color=vline_color, linewidth=2, linestyle=':',
        label='Median'
        )
    _ = g.ax_joint.axhline(
        y=df_stats['frag_mean'], xmin=0, xmax=1,
        color=vline_color, linewidth=3, linestyle='--',
        )
    _ = g.ax_joint.axhline(
        y=df_stats['frag_median'], xmin=0, xmax=1,
        color=vline_color, linewidth=3, linestyle=':',
        )

    g.fig.set_figwidth(20)
    g.fig.set_figheight(10)

    # Build legend for mean and median
    g.ax_joint.legend(
        loc='lower left',
        fontsize=24,
        markerscale=5,
        numpoints=20,
        frameon=False,
        ncol=2
        )

    # adjust layout, save, and close
    #plt.gcf().set_tight_layout(True)

    g.savefig(f'{op}_correlation_plot.png')
    plt.close()


def gather_data(ani):
    """Reads in the tsv ALLvsAll ANI file from MiGA/09.distances/03.ani

    Parameters
    ----------
    ani : str
        This is the miga-project.txt.gz file from the above directory

    Returns
    -------
    df : DataFrame
        DataFrame with two columns 'xs' and 'ys' containing all data
        This DataFrame also has a Comp column used to color code by
        category. For this option you need to adjust the A,B, n1 and n2
        parameters in this function to fit your delimiting/naming scheme.
    """

    data_dict = {'xs': [], 'ys': [], 'Comp': []}
    with open(ani, 'r') as f:
        header = f.readline()
        for l in f:
            X = l.rstrip().split('\t')
            n1 = X[1].split('_')[2]
            n2 = X[2].split('_')[2]

            Comp = f'{n1}-{n2}'

            ani = float(X[3])
            shared = float(X[5])
            total = float(X[6])
            ratio = shared / total

            data_dict['xs'].append(ani)
            data_dict['ys'].append(ratio)
            data_dict['Comp'].append(Comp)

    df = pd.DataFrame(data_dict)
    df = df[df['xs'] != 100.0]
    df = df[df['xs'] >= 90.0]

    return df

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-ani', '--ANI_tsv_file',
        help='Please specify the MiGA all vs all ani file!',
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
        '-op', '--output_prefix',
        help='What do you want to name the output file?',
        metavar='',
        type=str,
        required=True)
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('Running Script...')
    df = gather_data(
        args['ANI_tsv_file']
        )

    ANI_vs_shared_plot(
        df,
        args['organism_name'],
        args['output_prefix']
        )


if __name__ == "__main__":
    main()
