#/usr/bin/bash

EPA_Trees=$1
Output_Dir=$2

if [ ! -d ${Output_Dir} ]; then mkdir $Output_Dir; fi

for f in ${EPA_Trees}/*.jplace
  do
	x=`basename $f | cut -d. -f2`
	y=`echo $x | cut -d_ -f1-3`
	JPlace.to_iToL.rb -i ${f} -o ${Output_Dir}/$x -u $y
done

## To Run
## bash 00_PBS/07_Convert_EPA_Tree_JPlace.sh 06_EPA_Tree_Dir Output_Dir
## bash 00_PBS/07_Convert_EPA_Tree_JPlace.sh 06_Sruber_rpoB_EPA_Tree 07_Sruber_rpoB_iTol

