#/usr/bin/bash

Script_Path=/gpfs/pace2/project/bio-konstantinidis/rconrad6/Salinibacter_Assemblies/00_PBS/02b_Mash_sketch.pbs

inDir=$1
sffx=$2
outDir=$3
kmer=$4
hash_size=$5

query_array=(${inDir}/*.${sffx})
echo ${query_array[@]} > 02_${outDir}_query_array.txt
query_count=${#query_array[@]}
query_number=`expr $query_count - 1`

echo 'The Query count is' $query_count '. The Array count is' $query_number

if [ ! -d 00_log ]; then mkdir 00_log; fi
if [ ! -d ${outDir} ]; then mkdir ${outDir}; fi
if [ ! -d ${outDir}/Sketches ]; then mkdir ${outDir}/Sketches; fi

qsub -t 0-${query_number} -v outDir=${outDir},k=${kmer},s=${hash_size} ${Script_Path}


## Execute from within Blast Directory
##
## Example Command:
## bash 00_PBS/02a_Launch_Mash_sketch.sh 01_Sruber_Draft_Genomes fna 02_Mash_Analysis 31 1000000
## bash 00_PBS/02a_Launch_Mash_sketch.sh 01_Sruber_Draft_Genomes_Trim fna 02_Mash_Analysis 31 1000000
