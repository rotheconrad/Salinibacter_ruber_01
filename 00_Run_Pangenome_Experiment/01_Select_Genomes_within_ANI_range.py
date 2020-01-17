#!~/.conda/envs/rothon/bin/python

## USAGE :: python scriptname.py -gd Genomes_Directory -org Organisms -atl Lower_ANI_Threshold -atu Upper_ANI_Threshold
## Randomly Selects Genomes from a directory that are within an ANI_Threshold of eachother.
## Author :: Roth Conrad :: rotheconrad@gatech.edu :: https://github.com/rotheconrad
## Date Created :: Thursday, June 20, 2019
## Date Updated :: N/A first version

## python Scripts/Select_Genomes_above_ANI.py -n 5 -gd TEST_SAMPLES/ -org E_coli -atl 96 -atu 99

def choice_continue():
	yes = {'yes','y', 'ye', ''}
	no = {'no','n'}
	choice = input("Would you like to try again? [y/n]").lower()
	if choice in yes:
		return True
	elif choice in no:
		return False
	else:
		sys.stdout.write("Please respond with 'yes / y' or 'no / n'")

def if_gzip(new_dir, name):
	_ = subprocess.run(['gunzip', f'{new_dir}{name}'])
	name = '.'.join(name.split('.')[:-1])
	return name

def setup_the_bomb(e, n, gd, org):
	''' Creates new folder, dictionary, list of directories in genome directory, chooses g_prime, sets up the bomb '''

	# initialize a dictionary to keep a key of new files names with values of old file names
	file_name_key = {} 
	# get the list of files from the genome directory
	g_list = [g for g in os.listdir(gd) if os.path.isfile(f'{gd}{g}')]
	# make a new directory for this sampling
	new_dir = f'PGE_{e}_{org}/'
	_ = subprocess.run(['mkdir', new_dir])
	# Select the first genome of this sample, copy, unzip, and rename it
	g_prime = random.choice(g_list)	
	_ = g_list.remove(g_prime)
	_ = subprocess.run(['cp', f'{gd}{g_prime}', new_dir])
	if g_prime.split('.')[-1] == 'gz': g_prime = if_gzip(new_dir, g_prime)
	new_name = f'{org}_randomG_0001.fna'
	_ = subprocess.run(['mv', f'{new_dir}{g_prime}', f'{new_dir}{new_name}'])
	# Update the file name dictionary
	file_name_key[new_name] = g_prime

	return file_name_key, g_list, new_dir, new_name

def select_genomes(e, n, gd, org, atl, atu):
	'''
	Selects n genomes randomly from genome directory that are above the ANI threshold.
	Renames them iteratively and copies them to a n_organism_name_samples directory.
	Writes a key to match n sample to original genome file name.
	'''
	print('Selecting Random Genome 1 to use as ANI reference.')
	print('##################################################\n\n')
	file_name_key, g_list, new_dir, g_prime = setup_the_bomb(e, n, gd, org)
	genomes_avail = len(g_list)+1

	i, j = 1, 2
	while (j < n+1):
		if len(g_list) == 0:
			print('\n\n############# FAIL ####################################')
			print(f'Sampled all {genomes_avail} genomes in {gd} directory and only found {j} of the {n} matches requested between {atl}% and {atu}% ANI!\n\n')
			#select_genomes(n, gd, org, atl, atu) if choice_continue() else sys.exit()
			sys.exit()
		if i == 100000:
			print('\n\n############# FAIL ####################################')
			print(f'Sampled 100,000 Genomes and and only found {j} of the {n} matches requested between {atl}% and {atu}% ANI\n\n')
			sys.exit()
		i += 1
		print(f'Testing Random Genome {i} of {genomes_avail}')
		g_new = random.choice(g_list)
		print('g_list len:', len(g_list))
		_ = g_list.remove(g_new)
		_ = subprocess.run(['cp', f'{gd}{g_new}', new_dir])
		if g_new.split('.')[-1] == 'gz': g_new = if_gzip(new_dir, g_new)
		new_name = f'{org}_randomG_{j:04d}.fna'
		_ = subprocess.run(['mv', f'{new_dir}{g_new}', f'{new_dir}{new_name}'])
		fastANI = f"fastANI -r {new_dir}{g_prime} -q {new_dir}{new_name} -o {new_dir}temp.ani"
		out = subprocess.Popen(fastANI, shell=True, stdout=subprocess.PIPE)
		out = out.communicate()[0]
		if os.stat(f'{new_dir}temp.ani').st_size > 0:
			with open(f'{new_dir}temp.ani', 'r') as f:
				ani = float(f.readline().split('\t')[2])
				print(ani)

			if ani < atl:
				_ = subprocess.run(['rm', f'{new_dir}{new_name}'])
				print(f'FAIL: {ani}% ANI is below the {atl}% ANI threshold!')

			elif ani > atu:
				_ = subprocess.run(['rm', f'{new_dir}{new_name}'])
				print(f'FAIL: {ani}% ANI is above the {atu}% ANI threshold!')

			elif ani >= atl and ani < atu:
				file_name_key[new_name] = g_new
				print(f'PASS: {ani}% ANI is within the {atl}% to {atu}% threshold range!\nFound {j} genomes of the {n} requested.')
				j += 1

			else:
				print('FREAK OUT!!')

		else: subprocess.run(['rm', f'{new_dir}{new_name}'])

	_ = subprocess.run(['rm', f'{new_dir}temp.ani'])

	with open(f'{new_dir}00_Genome_Name_Key.tsv', 'w') as o:
		for k,v in file_name_key.items():
			o.write(f'{k}\t{v}\n')

	print('\n\n############# SUCCESS ####################################')
	print(f'{j-1} of {n} Genomes between {atl}% and {atu}% ANI sampled successfully!')
	print(f'Genome samples copied, Genome name key written to {new_dir}')
	print(f'Congratulations and enjoy your research.')


def main():

	# Configure Argument Parser
	parser = argparse.ArgumentParser(description='Randomly Selects Genomes from a directory that are within an ANI_Threshold of eachother.')
	parser.add_argument('-e', '--e_experiment', help='Experiment number', required=True)
	parser.add_argument('-n', '--n_genomes', help='How many genomes would you like to select?', type=int, required=True)
	parser.add_argument('-gd', '--genomes_directory', help='Please specify a directory of genomes to sample! (fasta or gzipped fasta)', required=True)
	parser.add_argument('-org', '--organism_name', help='What is the name of this organism?', required=True)
	parser.add_argument('-atl', '--ani_threshold_lower', help='What is the Lower ANI Threshold to use as the genome similarity cutoff?', type=float, required=True)
	parser.add_argument('-atu', '--ani_threshold_upper', help='What is the Upper ANI Threshold to use as the genome similarity cutoff?', type=float, required=True)
	args=vars(parser.parse_args())
	# Run this scripts main function
	print('\n\nRunning Script...\n\n')
	select_genomes(args['e_experiment'], args['n_genomes'], args['genomes_directory'], args['organism_name'], args['ani_threshold_lower'], args['ani_threshold_upper'])

if __name__ == "__main__":
	import argparse, sys, os, re, random, subprocess
	main()
