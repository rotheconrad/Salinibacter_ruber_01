#/usr/bin/bash

Script_Path=/gpfs/pace2/project/bio-konstantinidis/rconrad6/Salinibacter_Assemblies/00_PBS/05b_Blast_Check_Halorubrum.pbs

input_dir=$1
input_ext=$2
out_dir=$3
db=$4
pid=$5

file_array=(${input_dir}/${input_ext})
echo ${file_array[@]} > 05_Blast_query_list_${out_dir}.txt
file_count=${#file_array[@]}
file_number=`expr $file_count - 1`

if [ ! -d 00_log ]; then mkdir 00_log; fi
if [ ! -d $out_dir ]; then mkdir $out_dir; fi

echo Launching array of $file_count Blast jobs.

qsub -t 0-${file_number} -v out_dir=$out_dir,db=$db,pid=$pid ${Script_Path}

## Run Log:
## bash 00_PBS/05a_Launch_Blast_Check.sh 01_Sruber_Selected_Genomes *fna 05_Blast_Check_Blasts 01_Halorubrum_NCBI_Genomes/00_Halorubrum_DB.fna 90
