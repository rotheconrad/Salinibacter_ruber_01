#/usr/bin/bash

# Adjust script path
script=/gpfs/pace2/project/bio-konstantinidis/rconrad6/Salinibacter_Ramon/01_rpoB_pieTrees/00_PBS/03b_RefTree_RAxML.pbs

input_dir=$1
output_dir=$2

if [ ! -d ${output_dir} ]; then mkdir $output_dir; fi

input_array=(${input_dir}/*.aln)
input_count=${#input_array[@]}

echo Running RAxML on $input_count gene alignments.

qsub -t 0-$((input_count-1)) -v input_dir=$input_dir,output_dir=$output_dir $script

## To Run
## bash 00_PBS/03a_RefTree_Launch_RAxML.sh input_dir output_dir
## bash 00_PBS/03a_RefTree_Lauch_RAxML.sh 02_longSequence_Alignments 03_longSequence_RefTrees
