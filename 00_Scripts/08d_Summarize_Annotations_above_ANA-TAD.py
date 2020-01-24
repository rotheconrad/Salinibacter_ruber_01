#!/usr/bin/env python

'''

Plot Annotation/TAD Summary Results for file:
04_ClstrRepSeq_PanCat_ANATAD_Annotation.tsv

writes stacked barplots as .png files.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: Sunday, October 06, 2019
License :: GNU GPLv3
Copyright 2019 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
from collections import defaultdict
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pandas as pd

def get_summaries(annoTAD, TAD_thld):
    """ Read input table and generate counts """

    data = {
        'KEGG': defaultdict(list),
        'TrEMBL': defaultdict(list),
        'SwissProt': defaultdict(list)
        }

    df = pd.read_csv(annoTAD, sep='\t', index_col=0)

    PanCat = ['Core', 'Common', 'Rare', 'Specific']
    db = ['KEGG', 'TrEMBL', 'SwissProt']

    # gn is also the row names or index
    gn = [
        'hypothetical', 'uncharacterized', 'transferase', 'transposase',
        'synthase', 'integrase', 'transporter', 'recombination', 'synthesis',
        'transcriptase', 'reductase', 'kinase', 'regulator'
        ]

    colnames = [
                'Total_High', 'Total_Low', 'Core_High', 'Core_Low',
                'Common_High', 'Common_Low', 'Rare_High', 'Rare_Low',
                'Specific_High', 'Specifi_Low'
                ]

    for d in db: # for each database

        other = df

        for g in gn: # for each gene category in gn list calculation for total

             # Select gene category
            geneD = other[other[d].str.contains(g, case=False, regex=True)]
            # Select remaining genes
            other = other[~other[d].str.contains(g, case=False, regex=True)]

            geneD_high = geneD[geneD['ANA-TAD'] >= TAD_thld] # High TAD genes
            geneD_low = geneD[geneD['ANA-TAD'] < TAD_thld] # Low TAD genes

            data[d][g].append(float(len(geneD_high)))
            data[d][g].append(float(len(geneD_low)))

            for p in PanCat: # Calculations for each pangenome category

                p_high = geneD_high[geneD_high['PanCat'] == p]
                p_low = geneD_low[geneD_low['PanCat'] == p]

                data[d][g].append(float(len(p_high)))
                data[d][g].append(float(len(p_low)))


    ### Calculate "other" gene category ###
        otherD_high = other[other['ANA-TAD'] >= TAD_thld] # High TAD genes
        otherD_low = other[other['ANA-TAD'] < TAD_thld] # Low TAD genes

        data[d]['other'].append(float(len(otherD_high)))
        data[d]['other'].append(float(len(otherD_low)))

        for p in PanCat: # Calculations for each pangenome category

            p_high = otherD_high[otherD_high['PanCat'] == p]
            p_low = otherD_low[otherD_low['PanCat'] == p]

            data[d]['other'].append(float(len(p_high)))
            data[d]['other'].append(float(len(p_low)))


    ### Convert Dictionaries to DataFrames ####
    KEGG_df = pd.DataFrame.from_dict(data['KEGG'], orient='index', columns=colnames)
    TrEMBL_df = pd.DataFrame.from_dict(data['TrEMBL'], orient='index', columns=colnames)
    sProt_df = pd.DataFrame.from_dict(data['SwissProt'], orient='index', columns=colnames)

    return KEGG_df, TrEMBL_df, sProt_df

def Plot_Annotation_Summary(df, post, out):
    """ Plot Stacked Bar Charts by PanCats for each DB """

    colors = [
            '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c', '#8c510a',
            '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', '#01665e'
            ]

    normed_df = df.div(df.sum(axis=0), axis=1)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20,15))

    ax1 = df.T.plot.bar(stacked=True, ax=ax1, color=colors)
    ax2 = normed_df.T.plot.bar(stacked=True, ax=ax2, color=colors)

    # set the axis parameters / style
    for ax in fig.axes:
        ax.minorticks_on()
        ax.tick_params(axis='both', labelsize=18)
        ax.tick_params(axis='x', labelrotation=45)
        # set grid style
        ax.yaxis.grid(
            which="minor", color='#f0f0f0', linestyle='--', linewidth=2.5
            )
        ax.yaxis.grid(
            which="major", color='#d9d9d9', linestyle='--', linewidth=3
            )
        ax.set_axisbelow(True)

    handles1, labels1 = ax1.get_legend_handles_labels()
    ax1.legend(
        handles1[::-1],
        labels1[::-1],
        loc='center left',
        bbox_to_anchor=(0.98, 0.5),
        fancybox=True,
        shadow=True,
        fontsize=18
        )
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(
        handles2[::-1],
        labels2[::-1],
        loc='center left',
        bbox_to_anchor=(0.98, 0.5),
        fancybox=True,
        shadow=True,
        fontsize=18
        )

    # adjust layout, save, and close
    fig.set_tight_layout(True)
    plt.savefig(f'{out}_{post}.png')
    plt.close() 


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-a', '--ANATAD_annotation_file',
        help='The 04_ClstrRepSeq_PanCat_ANATAD_Annotation.tsv file.',
        metavar='',
        type=str,
        #required=True
        )
    parser.add_argument(
        '-t', '--TAD_Threshold',
        help='Threshold to select genes above ANA-TAD. (ex: 0.25)',
        metavar='',
        type=float,
        #required=True
        )
    parser.add_argument(
        '-o', '--out_file',
        help='What do you want to name the output file?',
        metavar='',
        type=str,
        #required=True
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('Running Script...')

    KEGG_df, TrEMBL_df, sProt_df = get_summaries(
                                        args['ANATAD_annotation_file'],
                                        args['TAD_Threshold']
                                        )

    Plot_Annotation_Summary(KEGG_df, 'KEGG', args['out_file'])
    Plot_Annotation_Summary(TrEMBL_df, 'TrEMBL', args['out_file'])
    Plot_Annotation_Summary(sProt_df, 'SwissProt', args['out_file'])


if __name__ == "__main__":
    main()
