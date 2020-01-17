#/usr/bin/bash

# Don't forget to make blast dbs for the genomes
# Confirm the use of MagicBlast makeblastdb over Blast+ makeblastdb
# which makeblastdb
# for f in *.fna; do makeblastdb -dbtype nucl -in $f -out $f -parse_seqids; done

Script_Path=/gpfs/pace2/project/bio-konstantinidis/rconrad6/Salinibacter_Assemblies/00_PBS/07b_MagicBlast.pbs

# Metagenome fasta directory. Metagenomic Reads are the queries
query_dir=$1
# Path to Group 1 reference genomes.
# The reference genomes we make the blast databases for.
# These are the subjects or databases.
db_dir=$2
# output directory
out_dir=$3

# We will use the Group 1 genomes to launch the job array.
db_array=(${db_dir}/*.fna)
db_count=${#db_array[@]}
db_number=`expr $db_count - 1`

if [ ! -d ${out_dir} ]; then mkdir ${out_dir}; fi
if [ ! -d 00_log ]; then mkdir 00_log; fi

for f in ${query_dir}/*fastq
  do
	echo Running Magic Blast on metagenome `basename $f` against $db_count genomes
	qsub -t 0-${db_number} -v query=${f},db_dir=${db_dir},out_dir=${out_dir} ${Script_Path}
done

## Run Log:
## bash 07a_Launch_MagicBlast.sh metagenome_dir Refgenome_dir output_dir
## bash 00_PBS/07a_Launch_MagicBlast.sh 00_Metagenome_ReadsTrimmed 01_Sruber_Group1_Genomes 07_MagicBlast_Metagenomes_to_Group1
