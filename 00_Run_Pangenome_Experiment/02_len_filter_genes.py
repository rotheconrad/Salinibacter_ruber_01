#!/usr/local/pacerepov1/python/2.7/bin/python

## USAGE :: python scriptname.py file.fasta prefix
## Reads fasta file and replaces sequence names with prefix_# counting from 1 to the end.
## Writes file.ref matching original names to new names.

import sys, subprocess

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

def Fasta_rename_sequences(infile, len_filter):

#    X = infile.split('/')
#    prefix = X[1].split('.')[0]
#    outfile = X[0] + prefix + '.rename'
    outfile = infile + '.lenfilter'

    with open(infile, 'r') as f, open(outfile, 'w') as o:
        i = 0
        c = 0
        for name, seq in read_fasta(f):
            c += 1
            if len(seq) > len_filter:
                o.write(f'{name}\n{seq}\n')
                i += 1

    _ = subprocess.run(['mv', outfile, infile])
#    print(f'Infile: {infile}')
#    print(f'Outfile: {outfile}')
#    print(f'Prefix: {prefix}')
    print(f'#### {infile} ################################################')
    print(f'Kept {i} of {c} predicted gene sequences')
    print(f'Sequences below {len_filter} nucleotide filter cutoff: {c-i}')
    print(f'##############################################################\n\n')

def main():
    infile = sys.argv[1]
    len_filter = int(sys.argv[2])
    Fasta_rename_sequences(infile, len_filter)
    
if __name__ == "__main__":
    main()

