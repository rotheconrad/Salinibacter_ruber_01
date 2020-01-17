#!/usr/bin/env python

'''Clean up Row and Column Names

This script just cleans up the column and row
names for the all vs all mash distance matrix.

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

import argparse

def Cleanup_RowCol_Names(infile, outfile):
    
    with open(infile, 'r') as f, open(outfile, 'w') as o:

        # Grab the column names from the data table
        colNames = f.readline().rstrip().split('\t')[1:]

        # Setup list of new column names. First column is Sample names
        newColNames = ['Samples',]

        # Read through list of colNames and split out the short unique id
        # adjust this to the naming scheme of the project
        for n in colNames:
            newName = n.split('/')[1].split('_')[1]
            newColNames.append(newName)

        o.write('\t'.join(newColNames) + '\n')

        # Read through each line and change the row name
        # This also needs to be adjusted to the naming scheme of the project
        for l in f:
            row = l.split('\t')
            newName = row[0].split('/')[1].split('_')[1]
            row[0] = newName
            o.write('\t'.join(row))
            

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_file',
        help='Please specify the Mash Distance Matrix!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--out_file',
        help='Please specify the output file name!',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    print('\nRunning Script...\n')
    Cleanup_RowCol_Names(
                    args['input_file'],
                    args['out_file']
                    )

    print('\nLooks like the row and column names were successfully renamed!\n\n')

if __name__ == "__main__":
    main()
