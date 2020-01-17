#/usr/bin/bash

Script_Path=/gpfs/pace2/project/bio-konstantinidis/rconrad6/Salinibacter_Assemblies/00_PBS/00b_Trimmomatic_Array.pbs

file_dir=$1
file_ext=$2
out_dir=$3
file_array=(${file_dir}/${file_ext})
file_count=${#file_array[@]}
file_number=`expr $file_count - 1`

if [ ! -d 00_log ]; then mkdir 00_log; fi
if [ ! -d $out_dir ]; then mkdir $out_dir; fi


echo Launching array of $file_count Trimmomatic jobs.

qsub -t 0-${file_number} -v file_dir=$file_dir,file_ext=$file_ext,out_dir=$out_dir ${Script_Path}

## Run Log:
## bash 00_PBS/00a_Launch_Trimmomatic_Array.sh 00_DataC *fastq 00_DataTrimmed
