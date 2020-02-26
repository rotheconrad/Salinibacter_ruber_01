#!/usr/bin/env python

'''Get Representative Fasta Sequence for CD-HIT Cluster

This script was written to retrieve the fasta sequence of selected gene
clusters from the CD-HIT (https://github.com/weizhongli/cdhit/wiki)
output files.

It takes the *.fnn.clstr file and the representative sequences *.fnn file
from CD-HIT(-EST) along with a text file containing 1 cluster number per
line and returns a fasta file of the representative sequences for the
user-specified gene clusters.

This script takes the following input parameters:

    * clstr - CD-HIT cluster file (str)
    * fasta - CD-HIT representative sequences fasta file (str)
    * clist - Text file containing 1 cluster number per line. (str)
    * out - name for the output fasta file(str)

This script returns the following files:

    * fasta file containing the representative sequence for the requested
      gene clusters.

This script requires the following packages:

    * argparse

This file can also be imported as a module and contains the follwing 
functions:

    * get_seqs_in_clstr - reads CD-HIT clstr file and retrieves seq names
    * get_seq_fasta - reads fasta file and gets sequences
    * main - the main function of the script

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: February 26th, 2020
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


def get_seq_fasta(fasta, seqs, out):
    ''' Read fasta file and retrieve and the matching sequences '''

    with open(fasta, 'r') as f, open(out, 'w') as o:
        for name, seq in read_fasta(f):
            name = name.split(' ')[0]
            if name in seqs:
                o.write(f"{name}\n{seq}\n")


def get_clstr_rep_seqs(clstr, clist):
    ''' reads CD-HIT clstr file and retrieves seq names '''

    seqs = {}
    cnumbs = []

    with open(clist, 'r') as f:
        for l in f: cnumbs.append(l.rstrip().split(', ')[0])

    with open(clstr, 'r') as f:
        for l in f:

            if l.startswith('>'):
                cluster = 'Cluster_' + l.rstrip().split(' ')[1]

            if not l.startswith('>') and cluster in cnumbs:
                X = l.rstrip().split(', ')[1].split('... ')
                gene = X[0]
                isrep = X[1]
                if isrep == '*': seqs[gene] = ''

    return seqs


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-c', '--cluster_file',
        help='Please specify the CD-HIT Cluster File!',
        metavar='',
        type=str,
        )
    parser.add_argument(
        '-f', '--fasta_file',
        help='Please specify the CD-HIT representative sequence fasta file!',
        metavar='',
        type=str,
        )
    parser.add_argument(
        '-n', '--cluster_number_list',
        help='Please specify the file with selected cluster numbers!',
        metavar='',
        type=str,
        )
    parser.add_argument(
        '-o', '--output_file',
        help='Please specify the output fasta file name',
        metavar='',
        type=str,
        )
    args=vars(parser.parse_args())

    clstr = args['cluster_file']
    clist = args['cluster_number_list']
    fasta = args['fasta_file']
    out = args['output_file']

    # Retrieve the list of sequence names for cnumb cluster
    print(f'Retrieving representative sequence names for cluster {clist}...')
    seqs = get_clstr_rep_seqs(clstr, clist)

    # Retrieve the fasta sequences for seq names in seqs. Write fasta output.
    print(f'Retrieving fasta sequences and writing to {out}')
    get_seq_fasta(fasta, seqs, out)


if __name__ == "__main__":
    main()
