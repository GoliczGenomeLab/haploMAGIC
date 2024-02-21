##!/bin/bash

#============================================================================================================================
#title: mendel.ME.rate.sh
#description: calculate ME rate per individual, marker and family
#author: jmontero
#date: 2023-03-20
#version: 1.0.0
#usage: bash mendel.ME.rate.sh "prefix1 prefix2 prefix3"
#notes: the input must contain the full prefix of the PED and MAP files (not only the population ID)
#	overcomes the plink contrain for which parents can only have one sex in the population
#	sorts by families the PED files organized by generations
#	calculating ME rate in simulated populations of known GE rate allows calculating GE rate in real data set from ME rate
#	knowing GE rate helps selecting the best haploMAGIC options by means of the simulation
#=============================================================================================================================

# Define global variables
pops=($1)
number_individuals=$(cat ${pops[0]}.ped | wc -l)

# Iterate over populations
for pop in "${pops[@]}"
do
	rm -f ${pop}_rfm.ped ${pop}_family
	awk '{$5="1"; print}' ${pop}.ped | tr ' ' '	' > ${pop}_newsex.ped # Create tmp PED with new sex column
	if [ ! -d PARENTS ]
	then
        	grep -v 'G0\|G1' ${pop}_newsex.ped | cut -f 3-4 | uniq > PARENTS # Create tmp file with parent column
	fi
	cat PARENTS | while read line # PARENTS contains the 3 and 4 th columns (father        mother)
	do
        	echo "$line" > tmp_family
	        family=$(cat tmp_family | tr '	' ':') # Create label for family
	        father=$(cut -f1 tmp_family) # Save father's name in var
        	mother=$(cut -f2 tmp_family) # Save mother's name in var
	        awk -v pat="$father" '$2 ~ pat { print $0 }' ${pop}_newsex.ped | tr ' ' '	' >> tmp_subset # Save the father's row in a tmp file
	        awk -v pat="$mother" '$2 ~ pat { print $0 }' ${pop}_newsex.ped | awk '$5="2" { print $0 }' | tr ' ' '	' >> tmp_subset # Save the mother's row in the tmp file and change the sex >
	        grep "$line" ${pop}_newsex.ped | grep "$mother" >> tmp_subset # Add the individuals of the family
	        generation=$(tail -n1 tmp_subset | cut -f1) # Extract the gen of the family
	        echo "$generation       $family" >> ${pop}_family # Fill up a file that links family with the children's generation
	        sed -E "s/^G[0-9]/"$family"/g" tmp_subset | sponge tmp_subset # Change the family label
	        cat tmp_subset >> ${pop}_rfm.ped # Add the tmp file to a reformatted file
		if [ ! -f ${pop}_rfm.map ]
		then
			cp ${pop}.map ${pop}_rfm.map
		fi
	        rm tmp*
	done
	rm -f PARENTS *newsex*
	plink1.9 --file ${pop}_rfm --mendel 1>/dev/null
	cut --complement -f 1-2 -d ' ' plink.lmendel | sed -E 's/ +/ /g' | sed -E 's/^ +//g' | sed 1d > ${pop}.marker.mendel
	cut --complement -f 1-2 -d ' ' plink.imendel | sed -E 's/ +/ /g' | sed -E 's/^ +//g' | sed 1d > ${pop}.individual.mendel
	cut --complement -f 1-2 -d ' ' plink.fmendel | sed -E 's/ +/ /g' | sed -E 's/^ +//g' | sed 1d > ${pop}.family.mendel
	cut --complement -f 1-2 -d ' ' plink.mendel | sed -E 's/ +/ /g' | sed -E 's/^ +//g' | sed 1d > ${pop}.info.mendel
	rm -f plink.*mendel *${pop}_rfm* ${pop}_family
done
rm -f plink.log
