#!/bin/bash

#==================================================================================
#title: 0.haplomagic.sh
#description: controls the execution of the pipeline
#author: jmontero
#date: 2023-03-20
#version: 1.0.0
#usage: bash 0.haplomagic.sh "pop1 pop2 pop3" "chr1 chr2 chr3"
#notes: Recommended options as default. Can be changed by rewriting the lines below
#==================================================================================

# Defining global variables
pops=($1)
chrs=($2)
min=2 # Options: any positive int (1/2/3/4...)
imp=None # Options: None/THonly/fullCor
cor=noCor # Options: noCor/fullCor/zeroOutFalseHom
thr=0 # Options: any positive int

# Recording execution time
runid=`date +%s`

# Executing pipeline
for pop in ${pops[@]} ;
do
	if [ -z $chrs ] ; then chrs=($(ls ${pop}_*.map | grep -oe '_[0-9].map' | sed 's/_//g' | sed 's/.map//g' | tr '\n' ' ')) ; fi # If no chrs provided, read all
	for chr in ${chrs[@]} ;
	do
		(chrstart=`date +%s` ;
		echo "pop=${pop}	chr=${chr}	state=1/5 phasing G0&G1" ; Rscript 1.phase.founders.R $pop $chr ;
		echo "pop=${pop}	chr=${chr}	state=2/5 phasing >G2" ; Rscript 2.phase.rest.R $pop $chr $min $imp $cor ;
		echo "pop=${pop}	chr=${chr}	state=3/5 haploblocking" ; Rscript 3.haploblock.R $pop $chr ;
		echo "pop=${pop}	chr=${chr}	state=4/5 detecting recombination" ; Rscript 4.detect.recombination.R $pop $chr $thr ;
		echo "pop=${pop}	chr=${chr}	state=5/5 calculating statistics" ; Rscript 5.statistics.R $pop $chr ;
		rm -f ${pop}_${chr}*{pedigree,genotype} ;
		chrend=`date +%s` ;
		chrruntime=$((chrend-chrstart)) ;
		echo "pop=${pop}	chr=${chr}	state=done (${chrruntime}s)") &>> haploMAGIC_${runid}_${pop}_${chr}.log ;
	done ;
done
