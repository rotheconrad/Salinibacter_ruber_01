#!~/.conda/envs/rothon/bin/python

## USAGE :: python scriptname.py -i infile.name -o outfile.name
## Insert short description of script purpose
## Author :: Roth Conrad :: rotheconrad@gatech.edu :: https://github.com/rotheconrad
## Date Created :: Insert date script created
## Date Updated :: N/A first version

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

def func1(f1, f2, outfile):
	
	d = {} # dictionary of sequence ids
	dd = {}


	with open(f1, 'r') as f:
		for name, seq in read_fasta(f):
			d[name] = ''

	with open(f2, 'r') as f, open(outfile, 'w') as o:
		for name, seq in read_fasta(f):
			if name in d:
				if seq[-1] == '*': seq = seq[:-1]
				o.write(name + '\n' + seq + '\n')

				if name in dd:
					dd[name] += 1
					#print(name)
				else:
					dd[name] = 0


	ds = {k:v for (k,v) in dd.items() if v > 0}
	print('Number of genes found more than once:', len(ds), '# This should be 0 or something is wrong!')

def main():

	# Configure Argument Parser
	import argparse
	parser = argparse.ArgumentParser(description='Grap sequences in file 1 from file 2.')
	parser.add_argument('-f1', '--file_one', help='Please specify the nucleotide sequence file!', required=True)
	parser.add_argument('-f2', '--file_two', help='Please specify the amino acid sequence  file!', required=True)
	parser.add_argument('-o', '--out_file', help='What do you want to name the output file?', required=True)
	args=vars(parser.parse_args())

	# Run this scripts main function
	print('Running Script...')
	func1(args['file_one'],args['file_two'], args['out_file'])

if __name__ == "__main__":
	main()
