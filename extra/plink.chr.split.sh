#!/bin/bash

#==========================================================================================
#title: plink.chr.split.sh
#description: split and recode ATGC PED and MAP files to 012 chromosome-based PED and MAP
#author: jmontero
#date: 2023-03-20
#version: 1.0.0
#usage: bash plink.chr.split.sh "pop1 pop2 pop3" "chr1 chr2 chr3"
#==========================================================================================

pops=($1)
chrs=($2)

for chr in ${chrs[@]} ;
do
	for pop in ${pops[@]} ;
	do
		plink1.9 --file $pop --chr $chr --recode 12 --out ${pop}_${chr} &>/dev/null
		rm -f ${pop}_${chr}.{nosex,log}
	done
done
