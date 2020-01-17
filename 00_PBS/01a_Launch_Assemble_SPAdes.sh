#/usr/bin/bash

switch=$4

if [ $switch -eq 0 ]
  then
Script_Path=/gpfs/pace2/project/bio-konstantinidis/rconrad6/Salinibacter_Assemblies/00_PBS/01b_Assemble_SPAdes.pbs

elif [ $switch -eq 1 ]
  then
Script_Path=/gpfs/pace2/project/bio-konstantinidis/rconrad6/Salinibacter_Assemblies/00_PBS/01c_Assemble_SPAdes.pbs

elif [ $switch -eq 2 ]
  then
Script_Path=/gpfs/pace2/project/bio-konstantinidis/rconrad6/Salinibacter_Assemblies/00_PBS/01d_Assemble_SPAdes.pbs

else
echo "Set switch=$4 to 0 for untrimmed reads, 1 for trimmomatic paired only reads, or 2 for trimmomatic paired/unpaired"
exit 1
fi

sample_list=$1
data_dir=$2
out_dir=$3
file_array=(`cat $sample_list`)
file_count=${#file_array[@]}
file_number=`expr $file_count - 1`

if [ ! -d 00_log ]; then mkdir 00_log; fi
if [ ! -d $out_dir ]; then mkdir $out_dir; fi


echo Launching array of $file_count SPAdes Assembly jobs.

qsub -t 0-${file_number} -v sample_list=$sample_list,data_dir=$data_dir,out_dir=$out_dir ${Script_Path}

## Run Log:
## bash 00_PBS/01a_Launch_Assemble_SPAdes.sh 01_Sample_List.txt 00_DataC 01_SPAdes_Assemblies 0
## bash 00_PBS/01a_Launch_Assemble_SPAdes.sh 01_Sample_List.txt 00_DataTrimmed 01_SPAdes_Assemblies_AdaptorTrimmed 1
## bash 00_PBS/01a_Launch_Assemble_SPAdes.sh 01_Sample_List.txt 00_DataTrimmed 01_SPAdes_Assemblies_Trim 2
