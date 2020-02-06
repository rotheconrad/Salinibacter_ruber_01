#/usr/bin/bash

# Adjust Path to Script
script=/gpfs/pace2/project/bio-konstantinidis/rconrad6/Salinibacter_Ramon/01_rpoB_pieTrees/00_PBS/05b_Align_ShortReads_MAFFT_reads.pbs

# Input Parameters
Short_Reads=$1
Ref_Aln=$2
Output_Dir=$3
sID=`basename $Ref_Aln | cut -d. -f1`

file_array=(${Short_Reads}/*.fasta) # adjust file extension to metagenomes
file_count=${#file_array[@]}

echo 'Adding MAFFt Short Read Alignments to '$file_count' long sequence alignments.'

if [ ! -d ${Output_Dir} ]; then mkdir ${Output_Dir}; fi
if [ ! -d 00_log/05_Align_ShortReads/ ]; then mkdir 00_log/05_Align_ShortReads/; fi

qsub -t 0-$((file_count-1)) -v Short_Reads=${Short_Reads},Ref_Aln=${Ref_Aln},Output_Dir=${Output_Dir},sID=${sID} ${script}

## Example Command:
## bash 00_PBS/05a_Align_ShortReads_Launch_MAFFT_reads.sh 04_Short_Reads_Dir 02_longSequence_alignment/sequence.aln Output_Dir
## bash 00_PBS/05a_Align_ShortReads_Launch_MAFFT_reads.sh 04_Sruber_rpoB_ShortReads 02_longSequence_Alignments/Srube_rpoB_Cluster_171_renamed.fasta.clustalo.aln 05_Sruber_rpoB_ShortRead_Alignments
## bash 00_PBS/05a_Align_ShortReads_Launch_MAFFT_reads.sh 04_Hqr02_rpoB_ShortReads 02_longSequence_Alignments/Hqr_rpoB_2.fasta.clustalo.aln 05_Hqr02_rpoB_ShortRead_Alignments
