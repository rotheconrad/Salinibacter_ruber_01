#!/usr/bin/env python

'''Combine Annotations Results.

This script was written for the Salinibacter ruber diversity project to
match the annotation pipeline results with results from the ANA-TAD
pipeline. In the end we have a tsv file with columns: 
Gene_Name, Gene_Cluster, PanCat, n/N, ANA-TAD, Annotations.

This tool takes the following input files:

    * 06_Annotations/04_ClstrRepSeq_Annotations_NoMatch.tsv
    * 08_Gene_TAD_Annotation_PanCat/00_Genes_Clusters_PanCats.tsv
    * 08_Gene_TAD_Annotation_PanCat/02_PanCats_ANATADs_byCluster.tsv

This script returns the following files:

    * tsv file with the combined results.

This script requires the following packages:

    * argparse 
    * collections.defaultdict

This file can also be imported as a module and contains the follwing 
functions:

    * gather_data - processes files into dictionaries.
    * combine_results - orchestrates the combining of results.
    * main - the main function of the script.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: Friday, August, 9th, 2019
License :: GNU GPLv3
Copyright 2019 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
from collections import defaultdict


def gather_ano(file):
    """ processes and filters ano data into a dict. returns dict."""

    # {sequence_name: [U-ID,TrEMBL,T-KO,U-ID,SwissProt,S-KO,KEGG,K-KO]}
    d = defaultdict(list)

    with open(file, 'r') as f:
        _ = f.readline()
        for l in f:
            X = l.rstrip().split('\t')
            db = X[0]
            gnm = X[1]
            uid = X[2]
            if uid == 'n/a': uid = 'No_Match'
            annotation = X[7]
            ko = X[9]
            if ko == 'n/a': ko = 'No_Match'
            if db == 'TrEMBL' or db == 'SwissProt':
                d[gnm].extend([uid, annotation, ko])
            elif db == 'KofamScan':
                d[gnm].extend([annotation, ko])

    return d


def gather_gcp(file):
    """ process data from lst file into a dict. returns dict. """

    d = {}

    with open(file, 'r') as f:
        _ = f.readline()
        for l in f:
            X = l.rstrip().split('\t')
            gnm = X[0]
            clstr = X[1]
            d[gnm] = clstr

    return d


def gather_tad(file):
    """ process data from the TAD file into dict. returns dict """

    d = {}
    tpoint = []

    with open(file, 'r') as f:
        header = f.readline().rstrip().split('\t')
        tpoints = '\t'.join(header[7:10])
        for l in f:
            X = l.rstrip().split('\t')
            clstr = X[0]
            partA = '\t'.join(X[:4])
            partB = '\t'.join(X[7:])
            data = f'{partA}\t{partB}'
            #data = f'{X[0]}\t{X[1]}\t{X[2]}\t{X[3]}\t{X[7]}\t{X[8]}\t{X[9]}\t{X[-1]}'
            d[clstr] = data

    return d, tpoints


def combine_results(tad, ano, gcp, out):
    """ preocesses files and outputs combined results to tsv """

    # Parse the input files
    # The metagenome timepoints are not in order from the previous files.
    # Label with tA, tB, and tC. User can specify order in the plot.
    # {Cluster_Name: [Cluster, PanCat, n/N, GeneLen, tA, tB, tC, ANA-TAD]}
    tad_dict, tpoints = gather_tad(tad)
    # {RepSeqName: Cluster_Name}
    gcp_dict = gather_gcp(gcp) # {RepSeqName: data}
    # {RepSeqName: [U-ID,TrEMBL,T-KO,U-ID,SwissProt,S-KO,KEGG,K-KO]}
    ano_dict = gather_ano(ano)

    # Set the output header
    header = (
        f'RepSeqName\tCluster\tPanCat\tn/N\tAvgLen\t{tpoints}\tANA-TAD\t'
        'U-ID\tTrEMBL\tT-KO\tKEGG\tK-KO\tU-ID\tSwissProt\tS-KO'
         )

    # Write new file
    with open(out, 'w') as o:

        o.write(header + '\n')

        for RepSeqName, annotation in ano_dict.items():
            clstr = gcp_dict[RepSeqName]
            anos = '\t'.join(annotation)
            newLine = f"{RepSeqName}\t{tad_dict[clstr]}\t{anos}"
            o.write(newLine + '\n')

    # Function and Script End.

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-t', '--ANA_TAD_file',
        help='The 02_PanCats_ANATADs_byCluster.tsv file',
        metavar='',
        type=str,
        #required=True
        )
    parser.add_argument(
        '-a', '--Annotation_file',
        help='The 04_ClstrRepSeq_Annotations_NoMatch.tsv file.',
        metavar='',
        type=str,
        #required=True
        )
    parser.add_argument(
        '-p', '--Genes_Clusters_PanCats',
        help='The 00_Genes_Clusters_PanCats.tsv file.!',
        metavar='',
        type=str,
        #required=True
        )
    parser.add_argument(
        '-o', '--out_file',
        help='What do you want to name the output .tsv file?',
        metavar='',
        type=str,
        #required=True
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('Running Script...')
    combine_results(
                        args['ANA_TAD_file'],
                        args['Annotation_file'],
                        args['Genes_Clusters_PanCats'],
                        args['out_file']
                        )


if __name__ == "__main__":
    main()
