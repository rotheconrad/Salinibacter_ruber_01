#/usr/bin/bash

script=/gpfs/pace2/project/bio-konstantinidis/rconrad6/Salinibacter_Ramon/01_rpoB_pieTrees/00_PBS/06b_EPA_Tree_RAxML.pbs

# Original
ShortRead_Alns=$1
Ref_Tree=$2
Output_Dir=$3
sID=`basename $Ref_Tree | cut -d. -f1`

if [ ! -d ${Output_Dir} ]; then mkdir $Output_Dir; fi
if [ ! -d 00_log/06_EPA_Tree ]; then mkdir 00_log/06_EPA_Tree; fi

input_array=(${ShortRead_Alns}/*.aln)
input_count=${#input_array[@]}

echo Running RAxML EPA on $input_count gene alignments.

qsub -t 0-$((input_count-1)) -v Ref_Tree=${Ref_Tree},ShortRead_Alns=${ShortRead_Alns},Output_Dir=${Output_Dir},sID=${sID} ${script}


## To Run
## bash 00_PBS/06a_EPA_Tree_Launch_RAxML.sh 05_ShortRead_Alignments_Dir Ref_Tree.aln Output_Dir
## bash 00_PBS/06a_EPA_Tree_Launch_RAxML.sh 05_Sruber_rpoB_ShortRead_Alignments 03_longSequence_RefTrees/RAxML_bestTree.Srube_rpoB_Cluster_171_renamed.fasta.clustalo.aln.tree 06_Sruber_rpoB_EPA_Tree
## bash 00_PBS/06a_EPA_Tree_Launch_RAxML.sh 05_Hqr01_rpoB_ShortRead_Alignments 03_longSequence_RefTrees/RAxML_bestTree.Hqr_rpoB_1.fasta.clustalo.aln.tree 06_Hqr01_rpoB_EPA_Tree
## bash 00_PBS/06a_EPA_Tree_Launch_RAxML.sh 05_Hqr02_rpoB_ShortRead_Alignments 03_longSequence_RefTrees/RAxML_bestTree.Hqr_rpoB_2.fasta.clustalo.aln.tree 06_Hqr02_rpoB_EPA_Tree
