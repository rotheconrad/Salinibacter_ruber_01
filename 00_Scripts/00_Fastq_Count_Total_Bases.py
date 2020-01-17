#!/usr/bin/env python

'''Calculate bases in fastq

Reads through directory of fastq files and
prints the total base pairs per file to 
stdout and the total base pairs of the project.

filename\tbase pairs sequenced
filename\tbase pairs sequenced
total\tbase pairs sequenced

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

def evaluate_data(fdir):

    f_list = [f for f in os.listdir(fdir) if os.path.isfile(f'{fdir}/{f}')]
    
    total = 0

    for i, file in enumerate(f_list):

        #print(f'{file}: {i+1:03}')
        line_count = 0
        subtotal = 0

        with open(f'{fdir}/{file}', 'r') as f:
            for l in f:
                line_count += 1
                if line_count%4 == 0:
                    slen = len(l.rstrip())
                    total += slen
                    subtotal += slen

        print(f'{file}\t{subtotal}')

    print(f'Total Bases:\t{total}')

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
    args=vars(parser.parse_args())

    fdir = args['file_directory']

    if fdir[-1] == '/': fdir = fdir[:-1]

    evaluate_data(fdir)

if __name__ == "__main__":
    main()
