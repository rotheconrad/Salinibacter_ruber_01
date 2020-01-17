#!/usr/bin/env python

'''Combine Annotations Results.

This script was written for the Salinibacter ruber diversity project to
match the annotation pipeline results with results from the TAD
by gene pipeline. In the end we have a list of the Gene Name / Gene 
Cluster with the Pangenome Category n/N, TAD value, and Annotation.

This tool takes the following input parameters:

    * The COMBINED_DBs_Annotation file from the annotation pipeline
    * The Gene_Cluster_Category_List file from the TAD by gene pipeline
    * One Avg_TAD_per_CLuster file from the TAD by gene pipeline
    * Minimum percent ID for an annotation.
    * Minimum Alignment / query length for an annotation.

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

def gather_lst(file):
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


def parse_kofam(X):
    """ Parses and filters lines for KofamScan DB """
    
    if X[7] == 'np.nan': annotation = 'Hypothetical Gene'
    else: annotation = X[7]

    return annotation


def parse_uniprots(X, pid, lng, uid):
    """ Parses and filters lines for UniProts DBs """

    aid = X[2]

    if aid == 'No_Match': annotation = 'Hypothetical Gene'

    else:
        pmatch = float(X[3])
        aratio = float(X[4]) / float(X[5])

        if pmatch >= pid and aratio >= lng:
            annotation = X[7]

        else:
            annotation = 'Hypothetical Gene'
            uid = 'No_Match'

    return annotation, uid


def gather_ano(file, pid, lng):
    """ processes and filters ano data into a dict. returns dict."""

    d = defaultdict(list)

    with open(file, 'r') as f:
        _ = f.readline()
        for l in f:
            if l.strip():
                X = l.rstrip().split('\t')
                db = X[0]
                gnm = X[1]
                uid = X[2]

                # Set the KO number or No_Match
                ko = X[9]
                if ko == 'np.nan' or ko == 'n/a': ko = 'No_Match'

                # Filter annotation results and return annotation
                if db == 'TrEMBL' or db == 'SwissProt':
                    annotation, uid = parse_uniprots(X, pid, lng, uid)
                    d[gnm].extend([uid, annotation, ko])

                elif db == 'KofamScan':
                    annotation = parse_kofam(X)
                    d[gnm].extend([annotation, ko])

                else: print('ERROR LINE:', l)

    return d


def gather_tad(file):
    """ process data from the TAD file into dict. returns dict """

    d = {}

    with open(file, 'r') as f:
        _ = f.readline()
        for l in f:
            X = l.rstrip().split('\t')
            clstr = X[0]
            d[clstr] = '\t'.join(X[1:])

    return d


def combine_results(tad, ano, lst, pid, lng, out):
    """ preocesses files and outputs combined results to tsv """

    # Parse the input files
    tad_dict = gather_tad(tad)
    ano_dict = gather_ano(ano, pid, lng)
    lst_dict = gather_lst(lst)

    # Set the output header
    header = (
        'Cluster\tName\tPanCat\tn/N\tAvgTAD\tAvgLen\tU-ID\tTrEMBL\t'
        'T-KO\tU-ID\tSwissProt\tS-KO\tKEGG\tK-KO\n'
         )

    # Write new file
    with open(out, 'w') as o:

        o.write(header)
        print(header)

        for k, v in ano_dict.items():
            clstr = lst_dict[k]
            anos = '\t'.join(v)
            newLine = f"{clstr}\t{k}\t{tad_dict[clstr]}\t{anos}\n"
            o.write(newLine)
            print(newLine)

    # Function and Script End.

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-tad', '--TAD_file',
        help='Avg_TAD_per_CLuster file from one metagenome.',
        metavar='',
        type=str,
        #required=True
        )
    parser.add_argument(
        '-ano', '--Annotation_file',
        help='Combined Annotation file.',
        metavar='',
        type=str,
        #required=True
        )
    parser.add_argument(
        '-lst', '--List_file',
        help='Gene_Cluster_Category_List file.!',
        metavar='',
        type=str,
        #required=True
        )
    parser.add_argument(
        '-pid', '--Percent_Match',
        help='Minimum Percent Identity for annotation. (eg 0.4)',
        metavar='',
        type=float,
        #required=True
        )
    parser.add_argument(
        '-lng', '--Match_Length',
        help='Minimum Alignment / Query length ratio. (eg 0.5)',
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
    combine_results(
                        args['TAD_file'],
                        args['Annotation_file'],
                        args['List_file'],
                        args['Percent_Match'],
                        args['Match_Length'],
                        args['out_file']
                        )


if __name__ == "__main__":
    main()
