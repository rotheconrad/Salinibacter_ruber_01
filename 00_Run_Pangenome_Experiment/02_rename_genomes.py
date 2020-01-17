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

def Fasta_rename_sequences(infile):

#    X = infile.split('/')
#    prefix = X[1].split('.')[0]
#    outfile = X[0] + prefix + '.rename'
    prefix = infile.split('.')[0]
    outfile = infile + '.rename'

    with open(infile, 'r') as f, open(outfile, 'w') as o:
        i = 1
        for name, seq in read_fasta(f):
            newName = '>%s_%d\n%s\n' % (prefix, i, seq)
            o.write(newName)
            i += 1

    _ = subprocess.run(['mv', outfile, infile])

    print(f'#### {infile} ################################################')
    print(f'Renamed {i-1} deflines.')
    print(f'##############################################################\n\n')

def main():
    infile = sys.argv[1]
    Fasta_rename_sequences(infile)
    
if __name__ == "__main__":
    main()

