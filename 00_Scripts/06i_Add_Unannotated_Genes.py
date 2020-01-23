#!/usr/bin/env python

'''Add UnAnnotated Genes to UniProt_Annotated.tsv files

Not all genes find a match to the UniProt DBs but its good to keep
track of them. This script adds no matches to the UniProt annotations.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: June 4th, 2019
License :: GNU GPLv3
Copyright 2019 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
from collections import defaultdict

def read_fasta(fp):
    ''' This function parses fasta files and yields name, seq'''
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


def Add_NoMatch_to_UniProt_Annotations(a, q, outfile):
    ''' This function adds unannotated genes to the output tsv'''

    # Build dictionary of {genes: {database: entry}}
    d = defaultdict(lambda: defaultdict())
    with open(a, 'r') as f:
            header = f.readline()
            for l in f:
                    X = l.split('\t')
                    database = X[0]
                    genename = X[1]
                    d[database][genename] = l

    # Read through fasta file of all representative genes
    # Write annotation for each gene for each database
    # if database annotation not available, write No_Match.
    with open(q, 'r') as f, open(outfile, 'w') as o:
        o.write(header)
        for name, seq in read_fasta(f):
            genename = name[1:].split(' ')[0]
            for database, genes in d.items():
                if genename in genes:
                    o.write(genes[genename])
                else:
                    o.write(
                        f'{database}\t{genename}\tn/a\tn/a\tn/a\t'
                        f'{len(seq)}\tn/a\tHypothetical Gene\tn/a\t'
                        f'n/a\tn/a\tn/a\tn/a\tn/a\tn/a\tn/a\n'
                        )


def main():

    # Configure Argument Parser
    import argparse
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-a', '--uniprot_annotated_tsv',
        help='Please specify the UniProt_Annotated.tsv file!',
        required=True,
        metavar='',
        type=str,
        )
    parser.add_argument(
        '-q', '--representative_protein_fasta',
        help='Please specify the representative protein fasta file!',
        required=True,
        metavar='',
        type=str,
        )
    parser.add_argument(
        '-o', '--out_file',
        help='What would you like to call the new Output file? (.tsv)',
        required=True,
        metavar='',
        type=str,
        )
    args=vars(parser.parse_args())

    # Run this scripts main function
    print('Running Script...')
    Add_NoMatch_to_UniProt_Annotations(
                            args['uniprot_annotated_tsv'],
                            args['representative_protein_fasta'],
                            args['out_file']
                            )

if __name__ == "__main__":
        main()

