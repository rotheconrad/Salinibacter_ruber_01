#!/usr/bin/env python

'''Filters 05 Blast Check Tables

Filters Blast matches for alignment length and percent identity.
Writes filtered output to input_file_name.filtered_best_hits.blst

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: January 14th, 2020
License :: GNU GPLv3
Copyright 2020 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse

def filter_blast(blast, alen, pID):

    outfile = blast.split('.')[0] + '.filtered_best_hits.blst'
    fails = 0 # counter for number of matches failing filters
    passes = 0 # counter for number of matches passing filters
    total = 0 # counter for total blast entries in file

    with open(blast, 'r') as file, open(outfile, 'w') as o:
        for line in file:
            total += 1
            X = line.rstrip().split('\t')
            alignment_length = int(X[3]) # read alignment length
            percent_ID = float(X[2]) # percent alignment of match

            if alignment_length >= alen and percent_ID >= pID:
                o.write(line)
                passes += 1
            else:
                fails += 1

    print('Total number of entries in blast file:', total)
    print('Number of reads failing the filters:', fails)
    print('Number of reads passing the filters:', passes)

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-b', '--blast_file',
        help='Please specify the tabular blast input file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-a', '--min_alignment_length',
        help='Minimum alignment length to keep in base pairs (ex: 700).',
        metavar='',
        type=int,
        required=True
        )
    parser.add_argument(
        '-p', '--min_percent_ID',
        help='Minimum percent identity to keep (ex: 90.0).',
        metavar='',
        type=float,
        required=True
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('Running Script...')
    filter_blast(
        args['blast_file'],
        args['min_alignment_length'],
        args['min_percent_ID']
        )


if __name__ == "__main__":
    main()

