#!/usr/bin/env python

'''Counts number of Base Pairs in a Fastq file.

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

def count_base_pairs(file):
  
    total = 0

    line_count = 0

    with open(file, 'r') as f:
        for l in f:
            line_count += 1
            if line_count%4 == 0:
                total += len(l.rstrip())

    return total

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--fastq_file',
        help='Please specify the fastq file!',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    infile = args['fastq_file']

    bps = count_base_pairs(infile)

    print(f'Total Base Pairs in {infile}:\t{bps}')

if __name__ == "__main__":
    main()