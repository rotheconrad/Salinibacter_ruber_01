#!/usr/bin/env python

'''Position(bp) vs Phred Score

Calculates the number of basepairs below the user defined quality
threshold across the dataset for each bp position of the reads.

Plots histogram of xaxis bp position and yaxis count of bp's below
threshold at that position across entire dataset.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: November 22nd, 2019
License :: GNU GPLv3
Copyright 2019 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse, os
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

def evaluate_data(fdir):

    f_list = [f for f in os.listdir(fdir) if os.path.isfile(f'{fdir}/{f}')]
    
    data = {
            'total': 0,
            'x<50': 0,
            '50<=x<75': 0,
            '75<=x<100': 0,
            '100<=x<150': 0,
            '150<=x': 0
            }

    for i, file in enumerate(f_list):

        print(f'{file}: {i+1:03}')
        line_count = 0

        with open(f'{fdir}/{file}', 'r') as f:
            for l in f:
                line_count += 1
                if line_count%4 == 0:
                    seq = l.rstrip()
                    seqlen = len(seq)
                    data['total'] += seqlen
                    if seqlen < 50:
                        data['x<50'] += seqlen
                    elif seqlen >= 50 and seqlen < 75: 
                        data['50<=x<75'] += seqlen
                    elif seqlen >= 75 and seqlen < 100:
                        data['75<=x<100'] += seqlen
                    elif seqlen >= 100 and seqlen < 150:
                        data['100<=x<150'] += seqlen
                    elif seqlen >= 150:
                        data['150<=x'] += seqlen

    return data


def plot_data(data, out):

    # Set Colors
    gridM = '#bdbdbd'
    gridm = '#d9d9d9'
    
    # prepare to plot data
    x = []
    y = []
    p = []

    x_order = [
                'total',
                'x<50',
                '50<=x<75',
                '75<=x<100',
                '100<=x<150',
                '150<=x'
                ]

    for i in x_order:
        x.append(i)
        y.append(data[i])
        p.append(f"{data[i]/data['total']*100:.2f}%")
        print(i, data[i])

    # Build the Plot
    fig, ax = plt.subplots(figsize=(15,10))

    # Plot titles and labels
    ax.set_title(f'Base pair counts by read lengths')
    ax.set_xlabel('x = Read Length')
    ax.set_ylabel(f'Total Base Pairs')

    # Setup up plot grids, spines and ticks
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

    # Plots
    ax.bar(x, y, width=0.8, align='center')

    for i, b in enumerate(ax.patches):
        width, height = b.get_width(), b.get_height()
        x, y = b.get_xy()
        ax.annotate(
            p[i],
            (x + width / 2, height),
            xytext=(0, 3),
            textcoords="offset points",
            ha='center', va='bottom',
            fontsize=24
            )

    # adjust layout, save, and close
    fig.set_tight_layout(True)
    plt.savefig(f'{out}_BPs_byLength.png')
    plt.close()


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-fdir', '--file_directory',
        help='Please specify the fastq file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-out', '--out_file',
        help='Please specify the output file prefix!',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    fdir = args['file_directory']
    out = args['out_file']

    if fdir[-1] == '/': fdir = fdir[:-1]

    data = evaluate_data(fdir)

    _ = plot_data(data, out)

if __name__ == "__main__":
    main()