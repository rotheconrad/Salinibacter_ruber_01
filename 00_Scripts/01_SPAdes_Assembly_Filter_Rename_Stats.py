#!/usr/bin/env python

'''Filters, Renames, and reports stats for SPAdes Assembly Contigs.

This script removes all contigs below the user specified value.
This script retains all contigs equal to and above the value.
This script removes all contigs with less than 2x coverage.
At the same time this script appends the user specified unique genome id
to the beginning of each contig name.
Finally, this script prints out a tab delimited line of stats:
ID\tTotal Contigs\tLength Fail\tCov Fail\tPass\tN50\tLength(bp)\t%GC

Pass is the number of contigs passing the filter and written to the
output file. The N50, Length(bp), and GC values are calculated for
the passing contigs.

This tool takes the following input parameters:

    * input file in fasta format. (reads, contigs, genes)

This script returns the following files:

    * output file of length filtered sequence in fasta format.

This script requires the following packages:

    * argparse

This file can also be imported as a module and contains the follwing 
functions:

    * read_fasta - parses fasta format to (name, seq)
    * print_lengths - prints length of each fasta entry to stdout.
    * length_filter - reads input, filters, writes output
    * main - the main function of the script

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: January 06, 2020
License :: GNU GPLv3
Copyright 2020 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse

def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

def print_lengths(infile):

    with open(infile, 'r') as f:
        for name, seq in read_fasta(f):
            print(f'{name}\t{len(seq)}')

def get_N50(length_list):
    half = sum(length_list) / 2
    descend_list = sorted(length_list, reverse=True)
    n = 0
    for i in descend_list:
        n += i
        if n >= half:
            return i

def get_GC(string):
    seq = string.upper()
    C = seq.count('C')
    G = seq.count('G')
    GC = (G + C) / len(seq) * 100
    return GC

def length_filter(infile, outfile, id, cutoff):

    total_count = 0
    pass_count = 0
    fail_count = 0
    cov_fail_count = 0
    cov_fails = []
    length_list = []
    seq_string = ''

    with open(infile, 'r') as f, open(outfile, 'w') as o:
        for name, seq in read_fasta(f):
            total_count += 1
            seq_len = len(seq)
            cov = float(name.split('_')[-1])

            if seq_len >= cutoff and cov >= 2:
                pass_count += 1
                length_list.append(seq_len)
                seq_string = seq_string + seq
                o.write(f'>{id}_{name[1:]}\n{seq}\n')

            elif seq_len < cutoff:
                fail_count += 1

            elif cov < 2:
                cov_fail_count += 1
                cov_fails.append(name[1:])

            else:
                print(f'Error with {name[1:]} length {seq_len}!')

    N50 = get_N50(length_list)
    total_length = len(seq_string)
    GC = get_GC(seq_string)

    stat_line = (
                    f'{id}\t'
                    f'{total_count}\t'
                    f'{fail_count}\t'
                    f'{cov_fail_count}\t'
                    f'{pass_count}\t'
                    f'{N50}\t'
                    f'{total_length}\t'
                    f'{GC:.2f}'
                )

    print(stat_line)

'''
    print(f'Sequences Longer or equal to {cutoff} bp:\t{pass_count}')
    print(f'Sequences Shorter than {cutoff} bp:\t{fail_count}\n')
    print(f'Sequences Long enough but less than 2x coverage: {len(cov_fails)}\n')
    for i in cov_fails:
        print(f'\t\t\t\t{i}')
    print('\n\n')
'''

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--in_file',
        help='Please specify an input file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--out_file',
        help='What do you want to name the output file?',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-g', '--unique_genome_id',
        help='Please specify a unique genome id!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-l', '--length_cutoff',
        help='What minimum length in bps do you want (int, ex: 1000)?',
        metavar='',
        type=int,
        required=True
        )
    parser.add_argument(
        '-p', '--print_lengths',
        help='Print lengths only (does not write file). (ex: -p 1)',
        metavar='',
        type=int,
        required=False,
        default=0
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
#    print('\nRunning Script...\n\n')
    if args['print_lengths'] == 0:
        length_filter(
                    args['in_file'],
                    args['out_file'],
                    args['unique_genome_id'],
                    args['length_cutoff']
                    )
    elif args['print_lengths'] == 1:
        print_lengths(
                    args['in_file']
                    )
    else:
        print('What is going on? Something is wrong.')


if __name__ == "__main__":
    main()

