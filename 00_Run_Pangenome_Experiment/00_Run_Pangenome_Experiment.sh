#/usr/bin/bash

## USAGE ####
## bash scriptname n_experiments n_genomes genome_directory Organism_name O_name ANI_lower ANI_upper Min_Gene_Len
## ie :: bash PGE.sh 100 102 Genome_Dir/ Salinibacter_ruber S_ruber 95 100 300
###################################################################################################
## This script launches all PBS Qsub jobs necessary for n number of experiments.
## Author :: Roth Conrad :: rotheconrad@gatech.edu :: https://github.com/rotheconrad
## Date Created :: Sunday, June 23, 2019
## Date Updated :: N/A first version
###################################################################################################
###################################################################################################

###############################################################################################
#################################################################################################
## Set Program Directory ####
###################################################################################################

PGE=/gpfs/pace2/project/bio-konstantinidis/rconrad6/Core-Pan_Genomes/00_Run_Pangenome_Experiment/

###############################################################################################
#################################################################################################
## Set variables to be input from command line in this order.
###################################################################################################

e=$1 # Number of experiments to run.
g=$2 # Number of genomes to sample per experiment.
gd=$3 # Directory with genomes in fasta or fasta.gz
org_name=$4 # Name of the organism like this: Escherichia_coli - used for plot titles
org=$5 # Organism abreviation like this E_coli - used for plot labels
atl=$6 # ANI Threshold Lower ( int or float ie 95 or 95.7 )
atu=$7 # ANI Threshold Upper ( int or float ie 100 or 98.9 )
len_filter=$8

log=PGE_log_${org}/
spath=/gpfs/pace2/project/bio-konstantinidis/rconrad6/Core-Pan_Genomes/00_Run_Pangenome_Experiment/
sd=PGE_0000_Summary/
if [ ! -d ${sd} ]; then mkdir ${sd}; fi
if [ ! -d ${sd}/A_${org}_Rarefaction_Summary/ ]; then mkdir ${sd}/A_${org}_Rarefaction_Summary/; fi
if [ ! -d ${sd}/B_${org}_Rarefaction_Results/ ]; then mkdir ${sd}/B_${org}_Rarefaction_Results/; fi
if [ ! -d ${sd}/C_${org}_Correlation_Summary/ ]; then mkdir ${sd}/C_${org}_Correlation_Summary/; fi
if [ ! -d ${sd}/D_${org}_Dereplicated_Genes/ ]; then mkdir ${sd}/D_${org}_Dereplicated_Genes/; fi

###################################################################################################
###################################################################################################

## This pipeline is split into separate PBS scripts because some steps depend on the completion 
## of the previous step and some steps are batch jobs that would overload the PBS server.

###################################################################################################
## Set directory for log output and echo variables
###################################################################################################

if [ ! -d ${log} ]; then mkdir ${log}; fi

printf "\n\nRunning Pangenome Experiments with the following settings:\n\n"

echo Number of Experiments to Run: $e
echo Number of Genomes to Sample: $g
echo Genome Directory to Sample from: $gd
echo Name of the Organism to Sample: $org_name
echo Abreviation of Organism Name: $org
echo Lower ANI Threshold value: $atl
echo Upper ANI Threshold value: $atu
echo Legnth Filter cutoff for gene predictions: $len_filter nucleotides
echo Directory for Log Output: $log
echo Path to Pangenome Experiment Scripts Directory: $spath

###################################################################################################
## Check the current state of the experiment
###################################################################################################

exp_total=$((e*g))
clstr_total=${e}
n_genomes=`echo PGE_0*_${org}/*fna | wc -w`
n_genes=`echo PGE_0*_${org}/03_FNN/*fnn | wc -w`
n_ani=`echo PGE_0*_${org}/00_ANI/*ani | wc -w`
n_clstr=`echo PGE_0*_${org}/*clstr | wc -w`
n_binary=`echo PGE_0*_${org}/*binary.tsv | wc -w`
n_correlation=`echo PGE_0*_${org}/*correlation_plot.pdf | wc -w`
n_rarefaction=`echo PGE_0*_${org}/*rarefaction_plot.pdf | wc -w`

# for tests to work when only running a single experiment.
if [ $n_genomes -eq 1 ] && [ ! -d PGE_0001_${org}/00_Genome_Name_Key.tsv ]; then n_genomes=0;fi
if [ $n_genes -eq 1 ] && [ ! -d PGE_0001_${org}/03_FNN/ ]; then n_genes=0; fi
if [ $n_ani -eq 1 ] && [ ! -d PGE_0001_${org}/00_ANI/ ]; then n_ani=0; fi
if [ $n_clstr -eq 1 ] && [ ! -s PGE_0001_${org}/*.clstr ]; then n_clstr=0; fi
if [ $n_clstr -eq ${e} ]; then n_genomes=${n_clstr}; n_genes=${n_clstr}; fi
if [ $n_clstr -lt ${e} ] && [ $n_clstr -ge 1 ]; then n_genomes=$((n_genomes+n_clstr)); n_genes=$((n_genes+n_clstr)); fi
if [ $n_binary -eq 1 ] && [ ! -s PGE_0001_${org}/*binary.tsv ]; then n_binary=0; fi
if [ $n_correlation -eq 1 ] && [ ! -s PGE_0001_${org}/*correlation_summary.tsv ]; then n_correlation=0; fi
if [ $n_rarefaction -eq 1 ] && [ ! -s PGE_0001_${org}/*rarefaction_results.tsv ]; then n_rarefaction=0; fi

printf '\n\nCurrent status of the Experiment:\n\n'

echo 'Total number of genomes this experiment will select:' $exp_total
echo 'Total genomes selected:' $n_genomes
echo 'Total genomes with genes predicted:' $n_genes
echo 'Total genomes with one vs all fastANI calculations:' $n_ani
echo 'Total experiments with gene clustering complete:' $n_clstr
echo 'Total experiments with binary gene cluster matrix:' $n_binary
echo 'Total experiments with ANI correlation results:' $n_correlation
echo 'Total experiments with pangenome rarefaction results:' $n_rarefaction

printf '\n\nStarting Experiment:\n\n'

                                                   ############################################ 
################################################################################################# 
################################################################################################### 
## Step 1 #### log file: PGE-Select_Genomes.out
## Iterate over experiment number and Randomly sample n genomes for each experiment.
## Places sampled genomes in PGE_${e_number}_${org}/ directory
###################################################################################################
##############

if [ ${n_genomes} -lt ${exp_total} ] && [ ${n_genes} -lt ${exp_total} ] && [ ${n_clstr} -lt ${clstr_total} ] && [ ${n_ani} -ne ${exp_total} ]
  then
	echo '01 - Selecting Genomes for each experiment.'
	s1=$(qsub -t 1-${e} -v spath=${spath},log=${log},n=${g},gd=${gd},org=${org},atl=${atl},atu=${atu} ${PGE}01_Select_Genomes.pbs)
	step1="-W depend=afterokarray:`echo $s1 | cut -d[ -f1`[]"
#	step1=`echo $s1 | cut -d[ -f1`[]
	echo $step1
	sleep 3
  else
	echo '01 - Looks like the genomes have all been selected.'
#	s1=$(qsub -t 1-1 -v log=${log} ${PGE}01_Select_Genomes_Completed.pbs)
#	step1=`echo $s1 | cut -d[ -f1`[]"
	step1=""
	echo $step1
	sleep 3
fi
                                                    ########################################### 
################################################################################################# 
################################################################################################### 
## Step 2 #### log file: PGE-predictGenes_launch.out
## Once step 1 is finished. Run gene predictions for all genomes within each experiment.
## Also runs all by all fastANI Calculations.
## Takes about ~2 hours for 100 genomes. Adjust walltime as appropriate.
## Updated pipeline to use GNU paralell. This Step now includes FastANI from Step 3.
## Step 3 is no longer needed.
###################################################################################################
##############

if [ ${n_genes} -lt ${exp_total} ] && [ ${n_clstr} -lt ${clstr_total} ] || [ ${n_ani} -lt ${exp_total} ]
  then
	echo '02 - Predicting the genes for each experiment.'
#	s2=$(qsub -W depend=afterokarray:${step1} -t 1-${e}%2 -v spath=${spath},log=${log},g=${g},org=${org},len_filter=${len_filter} ${PGE}02_predictGenes_launch.pbs)
	s2=$(qsub ${step1} -t 1-${e} -v spath=${spath},log=${log},g=${g},org=${org},len_filter=${len_filter} ${PGE}02_predictGenes_launch.pbs)
	step2="-W depend=afterokarray:`echo $s2 | cut -d[ -f1`[]"
	echo $step2
	sleep 3
  else
	echo '02 - Looks like the genes have all been predicted.'
#	s2=$(qsub -t 1-1 -v log=${log} ${PGE}02_predictGenes_Completed.pbs)
#	step2=`echo $s2 | cut -d[ -f1`[]
	step2=""
	echo $step2
#	sleep 3
fi

                                                    ###########################################
#################################################################################################
###################################################################################################
## Step 3 #### log file: PGE-fastANI_launch.out
## Once step 2 is finished. Run all by all fastANI on all genomes within each experiment.
###################################################################################################
##############

#if [ ${n_ani} -lt ${exp_total} ] && [ ${n_clstr} -lt ${clstr_total} ] && [ ${n_ani} -ne ${exp_total} ]
#  then
#	echo '03 - Calculating all vs all fastANI values for each experiment.'
#	s3=$(qsub -W depend=afterokarray:${step2} -t 1-${e}%2 -v spath=${spath},log=${log},org=${org} ${PGE}03_fastANI_launch.pbs)
#	step3=`echo $s3 | cut -d[ -f1`[]
#	echo $step3
#	sleep 3       
#  else
#	echo '03 - Looks like the fastANI values have all been calculated.'
#	s3=$(qsub -t 1-1 -v log=${log} ${PGE}03_fastANI_Completed.pbs)
#	step3=`echo $s3 | cut -d[ -f1`[]
#	echo $step3
#	sleep 3
#fi
                                             ########################################### 
################################################################################################# 
################################################################################################### 
## Step 4 #### log file: PGE-cluster_genes
## Once step 3 is finished.  Perfom gene clustering.
###################################################################################################
##############

if [ ${n_clstr} -lt ${clstr_total} ]
  then
	echo '04 - Clustering genes for each experiment.'
#	s4=$(qsub -W depend=afterokarray:${step3} -t 1-${e} -v log=${log},org=${org} ${PGE}04_cluster_genes.pbs)
	s4=$(qsub ${step2} -t 1-${e} -v log=${log},org=${org} ${PGE}04_cluster_genes.pbs)
#	s4=$(qsub -t 1-${e} -v log=${log},org=${org} ${PGE}04_cluster_genes.pbs)
#	step4=`echo $s4 | cut -d[ -f1`[]
	step4="-W depend=afterokarray:`echo $s4 | cut -d[ -f1`[]"
	echo $step4
	sleep 3
  else
	echo '04 - Looks like the genes have all been clustered.'
#	s4=$(qsub -t 1-1 -v log=${log} ${PGE}04_cluster_genes_Completed.pbs)
#	step4=`echo $s4 | cut -d[ -f1`[]
	step4=""
	echo $step4
	sleep 3
fi
                                                    ########################################### 
################################################################################################# 
################################################################################################### 
## Step 5 #### log file: PGE-cluster_analysis.out
## Once step 4 is finished.  Perfom gene clustering analysis.
###################################################################################################
##############

if [ ${n_binary} -lt ${clstr_total} ] || [ ${n_correlation} -lt ${clstr_total} ] || [ ${n_rarefaction} -lt ${clstr_total} ]
  then
	echo '05 - Analyzing clusters for each experiment.'
#	s5=$(qsub -W depend=afterokarray:${step4} -t 1-${e} -v log=${log},name=${org_name},prm=${g},spath=${spath},org=${org} ${PGE}05_cluster_analysis.pbs)
	s5=$(qsub ${step4} -t 1-${e} -v log=${log},name=${org_name},prm=${g},spath=${spath},org=${org} ${PGE}05_cluster_analysis.pbs)
	step5=`echo $s5 | cut -d[ -f1`[]
	echo $step5
	sleep 3
  else
	echo '05 - Looks like the clusters have all been analyzed.'
	sleep 3
fi
                                                    ########################################### 
################################################################################################# 
###################################################################################################

############################################ FINISH ###############################################

printf '\n\n ## Script END. The experiment should be running or it is complete.## \n\n'

###################################################################################################
#################################################################################################
                                                    ###########################################
