#!/usr/bin/env python

'''

Build Plots to compare annotation results for 2 difference organisms

writes stacked barplots as .png files.
writes other gene category annotations to tsv files for each database.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: January 23rd, 2020
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
import pandas as pd

def get_summaries(f1, f2, out):
    """ Read input table and generate counts """

    data = defaultdict(list)

    s1 = f1.split('_')[1]
    s2 = f2.split('_')[1]
    df1 = pd.read_csv(f1, sep='\t', index_col=0)
    df2 = pd.read_csv(f2, sep='\t', index_col=0)

    db = ['KEGG', 'TrEMBL', 'SwissProt']

    # gn is also the row names or index
    # This is a list of key words to look for in the gene annotations
    gn = [
        'hypothetical', 'uncharacterized', 'transferase', 'transposase',
        'synthase', 'integrase', 'transport', 'recombination', 'synthesis',
        'transcriptase', 'transcription', 'reductase', 'kinase', 'regulator',
        'repair', 'inner membrane', 'outer membrane', 'membrane', 'replication',
        'synthetase', 'recombinase', 'helicase', 'endonuclease', 'ribonuclease',
        'nuclease', 'isomerase', 'secretion', 'chaperone', 'export', 'capsid',
        'hydroxylase', 'interferase', 'CRISPR', 'ribosomal', 'ribosome',
        'phosphatase', 'dehydrogenase', 'hydrogenase', 'lipoprotein',
        'protease', 'phosphodiesterase', 'rhodopsin', 'ATP-binding',
        'GTP-binding', 'DNA-binding', 'plasmid', 'phage', 'cytoplasmic',
        'phosphorylase', 'monooxygenase', 'polymerase', 'esterase', 'epimerase',
        'desulfurase','hydrolase', 'antiporter', 'heat shock|hsp',
        'dehydratase', 'transaldolase', 'peptidase', 'oxidase'
        ]

    # This is the 3 databases x the 2 organisms
    colnames = [
                f'{s1}_KEGG', f'{s2}_KEGG', f'{s1}_TrEMBL',
                f'{s2}_TrEMBL', f'{s1}_SwissProt', f'{s2}_SwissProt',
                ]

    Master_other1 = df1
    Master_other2 = df2
    total1 = df1.shape[0]
    total2 = df2.shape[0]

    for d in db: # for each database

        # the other category is built by subtracting results of each
        # search from the initial databases. Everything remaining after
        # searching all of the gn list is the other category.
        other1 = df1
        other2 = df2

        print(f'\n\n{d}')

        for g in gn: # for each gene category in gn list calculation for total

            # Select gene category
            select1 = other1[d].str.contains(g, case=False, regex=True)
            select2 = other2[d].str.contains(g, case=False, regex=True)
            Mselect1 = Master_other1[d].str.contains(g, case=False, regex=True)
            Mselect2 = Master_other2[d].str.contains(g, case=False, regex=True)

            geneD1 = other1[select1]
            geneD2 = other2[select2]
            # Select remaining genes
            other1 = other1[~select1]
            other2 = other2[~select2]
            Master_other1 = Master_other1[~Mselect1]
            Master_other2 = Master_other2[~Mselect2]

            len1 = geneD1.shape[0]
            len2 = geneD2.shape[0]
            data[g].append(len1)
            data[g].append(len2)

            print(
                g, len1, f'{len1/total1*100:.2f}',
                len2, f'{len2/total2*100:.2f}'
                )


        ### Calculate "other" gene category ###
        data['other'].append(float(len(other1)))
        data['other'].append(float(len(other2)))

        ### Write other category for each db to file ###
        other1.to_csv(f'{out}_{s1}_{d}.tsv', sep='\t')
        other2.to_csv(f'{out}_{s2}_{d}.tsv', sep='\t')

    ### Write master other category to file ###
    Master_other1.to_csv(f'{out}_{s1}_Master.tsv', sep='\t')
    Master_other2.to_csv(f'{out}_{s2}_Master.tsv', sep='\t')

    ### Convert Dictionaries to DataFrames ####
    df = pd.DataFrame.from_dict(data, orient='index', columns=colnames)

    return df


def Plot_Annotation_Summary(df, out):
    """ Plot Stacked Bar Charts by PanCats for each DB """

    colors = [
        #'#a6cee3','#b2df8a','#fb9a99','#8c510a','#ff7f00','#6a3d9a','#b15928',
        #'#1f78b4','#33a02c','#e31a1c','#fdbf6f','#cab2d6','#ffff99','#01665e'
        '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c', '#8c510a',
        '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', '#01665e'
            ]

    normed_df = df.div(df.sum(axis=0), axis=1)

    fig, ax = plt.subplots(figsize=(20,12))

    #ax1 = df.T.plot.bar(stacked=True, ax=ax1, color=colors)
    ax = normed_df.T.plot.bar(stacked=True, ax=ax, color=colors)

    # set the axis parameters / style
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

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(
        handles[::-1],
        labels[::-1],
        ncol=3,
        loc='center left',
        bbox_to_anchor=(0.98, 0.5),
        fancybox=True,
        shadow=True,
        fontsize=18
        )

    # adjust layout, save, and close
    fig.set_tight_layout(True)
    plt.savefig(out)
    plt.close() 


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-f1', '--annotation_file_1',
        help='The 1st 04_*_Transformed_Annotations.tsv file.',
        metavar='',
        type=str,
        #required=True
        )
    parser.add_argument(
        '-f2', '--annotation_file_2',
        help='The 2nd 04_*_Transformed_Annotations.tsv file.',
        metavar='',
        type=str,
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

    df = get_summaries(
                args['annotation_file_1'],
                args['annotation_file_2'],
                args['out_file']
                )

    Plot_Annotation_Summary(df, args['out_file'])


if __name__ == "__main__":
    main()
