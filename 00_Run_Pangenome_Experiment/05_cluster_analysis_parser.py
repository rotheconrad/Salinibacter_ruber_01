#!~/.conda/envs/rothon/bin/python

## USAGE :: python scriptname.py -clstr cd-hit_output.clstr -op output_prefix
## This script parses a CD-HIT clstr file into a binary presence/absence matrix of Gene Clusters per Genome.
## Author :: Roth Conrad :: rotheconrad@gatech.edu :: https://github.com/rotheconrad
## Date Created :: Friday, June 07, 2019
## Date Updated :: N/A first version

def _binary(ClusterIDs,GenomeIDs):
	'''
	Turns ClusterIDs and GenomeIDs into the binary matrix
	'''

	print('Generating Binary Matrix ...')
	d3 = defaultdict(list)
	for k,v in GenomeIDs.items(): # read through genomeID's
		for i in ClusterIDs:

			if i in v: # if genomeID has the gene cluster add a 1.
				d3[k].append(1)

			else: # if genomeID does not have the gene cluster add a 0.
				d3[k].append(0)

	d3['Cluster'] = ClusterIDs

	return d3

	
def Generate_Binary_from_CD_HIT_Clstr(clstr, op):
	'''
	Creates the binary presence/absence matrix of Gene clusters per genome.
	'''

	# Initialize variables.
	ClusterIDs = []
	GenomeIDs = defaultdict(lambda: defaultdict(int))

	clusterID = None

	# Read through clstr file and populate initialized variables.
	print('Parsing CD-HIT Clstr File ...')
	with open(clstr, 'r') as c:
		for l in c:
			if l[0] == '>': # This denotes the beginning of a new cluster.
				if clusterID: # When cluster is set, populate the dictionary with current cluster ID.
					ClusterIDs.append(clusterID)
					clusterID = None # reset cluster.

				clusterID = l[1:8] + '_' + l[9:].rstrip() # name new cluster.

			elif l[0].isdigit(): # This denotes entry within a cluster.
				X = l.rstrip().split(' ') # Split the line by space.

				## THIS IS WHAT NEEDS TO BE ADJUSTED TO YOUR FILE NAMES ##########################
				genomeID = X[1].split('_')[3] # Select the genome ID. ############################
				## THIS IS WHAT NEEDS TO BE ADJUSTED TO YOUR FILE NAMES ##########################

				GenomeIDs[genomeID][clusterID] += 1 # Keep track of which clusters each genome is in and how many times.

			else: # There should be nothing else. Print fun warning statement if there is.
				print('########### !!! FREAK OUT !!! Something is wrong. Stop and assess !!!  #############')

		# need to populate d1 with stats for the last cluster.
		ClusterIDs.append(clusterID)

	d3 = _binary(ClusterIDs, GenomeIDs)

	# Save the Parsed Cluster data
	binary = pd.DataFrame(d3)
	n2 = op + '.tsv'
	binary.to_csv(n2, sep='\t')

	return binary, n2

def main():

	# Configure Argument Parser
	parser = argparse.ArgumentParser(description='This script parses a CD-HIT clstr file into a binary presence/absence matrix of Gene Clusters per Genome.')
	parser.add_argument('-clstr', '--cdhit_clstr_file', help='Please specify a CD-HIT clstr file to us for input!', required=True)
	parser.add_argument('-op', '--output_prefix', help='What prefix would you like to use for the output files?', required=True)
	args=vars(parser.parse_args())

	# Run this scripts main function
	print('Running Script...')
	binary, n2 = Generate_Binary_from_CD_HIT_Clstr(args['cdhit_clstr_file'], args['output_prefix'])
	print(f'Success!! Binary matrix written to {n2} file.')

if __name__ == "__main__":
	import argparse
	from collections import defaultdict
	import pandas as pd
	main()
