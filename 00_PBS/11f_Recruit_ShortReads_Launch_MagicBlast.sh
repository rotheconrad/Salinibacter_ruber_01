#/usr/bin/bash

# Don't forget to make blast dbs for the genomes
# makeblastdb -dbtype nucl -in  -out  -parse_seqids
# makeblastdb -dbtype nucl -in All_Drafts_E1_M15_S4_fltrd.blast.fasta -out All_Drafts_E1_M15_S4_fltrd.blast.fasta -parse_seqids

# Adjust Path to Scripts
pbsscript=/gpfs/pace2/project/bio-konstantinidis/rconrad6/Salinibacter_Ramon/01_rpoB_pieTrees/00_PBS/04b_Recruit_ShortReads_MagicBlast.pbs
filter=/gpfs/pace2/project/bio-konstantinidis/rconrad6/Salinibacter_Ramon/01_rpoB_pieTrees/00_Scripts/00_pieTree_MagicBlast_filter.py
collect=/gpfs/pace2/project/bio-konstantinidis/rconrad6/Salinibacter_Ramon/01_rpoB_pieTrees/00_Scripts/01_Get_FastaReads_Matching_Filtered_Blast.py

# Input Parameters
Query_Path=$1
Ref_Fasta=$2
sID=`basename $Ref_Fasta | cut -d. -f1`
Output_Dir=$3

# Build MagicBlast Database
echo Building Magic Blast database for $Ref_Fasta
makeblastdb -dbtype nucl -in $Ref_Fasta -out $Ref_Fasta -parse_seqids


query_array=(${Query_Path}/*.fa) # adjust file extension to metagenomes
query_count=${#query_array[@]}

echo 'Running Magic Blast for '$query_count' metagenomes.'


if [ ! -d ${Output_Dir} ]; then mkdir ${Output_Dir}; fi
if [ ! -d 00_log/04_Recruit_ShortReads/ ]; then mkdir 00_log/04_Recruit_ShortReads/; fi

qsub -t 0-$((query_count-1)) -v Query_Path=${Query_Path},Ref_Fasta=${Ref_Fasta},sID=${sID},Output_Dir=${Output_Dir},filter=${filter},collect=${collect} ${pbsscript}

## Example Command:
## bash 00_PBS/04a_Recruit_ShortReads_Launch_MagicBlast.sh metagenome_dir ref_fasta output_dir
## bash 00_PBS/04a_Recruit_ShortReads_Launch_MagicBlast.sh ../Saltern-metagenomes 01_Sruber_rpoB_sequence/Srube_rpoB_Cluster_171_renamed.fasta 04_Sruber_rpoB_ShortReads
