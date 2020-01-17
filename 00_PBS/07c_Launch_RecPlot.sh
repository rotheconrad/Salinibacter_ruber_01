#/usr/bin/bash

Script_Path=/gpfs/pace2/project/bio-konstantinidis/rconrad6/Core-Pan_Genomes/COVERAGE/PBS/02b_RecPlot.pbs

blast_dir=$1
genome_dir=$2
thrshld=$3


if [ ! -d ${blast_dir}Tab_Blasts/ ]; then mkdir ${blast_dir}Tab_Blasts; fi

if [ ! -d ${blast_dir}RecPlotPDF/ ]; then mkdir ${blast_dir}RecPlotPDF; fi

if [ ! -d ${blast_dir}RecPlotCatSBJ/ ]; then mkdir ${blast_dir}RecPlotCatSBJ; fi

if [ ! -d ${blast_dir}RecPlotRDATA/ ]; then mkdir ${blast_dir}RecPlotRDATA; fi

if [ ! -d ${blast_dir}log/ ]; then mkdir ${blast_dir}log; fi

query_array=(${blast_dir}*.blast)
echo ${query_array[@]} > ${blast_dir}log/RecPlot_temp_query_array.txt
query_count=${#query_array[@]}
query_number=`expr $query_count - 1`

echo 'The Query count is' $query_count '. The Array count is' $query_number

qsub -t 0-${query_number} -v blast_dir=$1,genome_dir=$2,thrshld=$3 ${Script_Path}


##############################################
# genome_dir - directory with reference genomes
# thrshld - sequence discrete population line: int or float ex: 95 or 97.5
# 95 is a good starting point.


## Execute from outside MiGA_project directory
## bash PBS/02a_Launch_RecPlot.sh Blast_Results/ genomes/ 95
## bash PBS/02a_Launch_RecPlot.sh MegaBlast_Results/ genomes/ 95
## bash PBS/02a_Launch_RecPlot.sh BLAT_Results/ genomes/ 95
## bash PBS/02a_Launch_RecPlot.sh BLAT-fine_Results/ genomes/ 95
## bash PBS/02a_Launch_RecPlot.sh Bowtie2-Blast_Results/ genomes/ 95
