#!/usr/bin/env Rscript

#===============================================================================================================================================
#title: simulate_population.R
#description: applies the AlphaSimR functions to generate a simulated population from input PED and MAP files with a given genotyping error rate
#author: jmontero
#date: 2023-03-20
#version: 1.0.0
#usage: Rscript simulate_population.R pop1 chr1 GErate
#===============================================================================================================================================


# Loading packages
suppressMessages(suppressWarnings(if(!require("hash")){ install.packages("hash") }))
suppressMessages(suppressWarnings(if(!require("AlphaSimR")){ install.packages("AlphaSimR") }))
suppressMessages(suppressWarnings(if(!require("readr")){ install.packages("readr") }))
suppressMessages(suppressWarnings(if(!require("stringr")){ install.packages("readr") }))

args = commandArgs(trailingOnly=TRUE)

# Define global parameters
popID=args[1]
chrID=args[2]
genErrorRate = as.numeric(args[3])/100
simID=paste0("simpop_", as.character(genErrorRate*100), "_",  popID, "_", chrID)

# Input pedigree file
pedDF = read.table(paste0(popID, "_", chrID, ".ped"), header = FALSE, sep = " ", colClasses = "character")
mapDF = read.table(paste0(popID, "_", chrID, ".map"), header = FALSE, sep = "\t")
realGenotypes = pedDF[, -c(1, 2, 3, 4, 5, 6)]
pedigreeInitial = pedDF[, c(1, 2, 3, 4)]
nFounders=sum(pedigreeInitial[,1] == 'G0')
firstG2 = which(pedigreeInitial[,1] == 'G2')[1]
firstG3 = which(pedigreeInitial[,1] == 'G3')[1]
firstG4 = which(pedigreeInitial[,1] == 'G4')[1]
nSnp = nrow(mapDF)

# Create inbred founder lines using the Markovian Coalescent Simulator
founderPop <- 
  runMacs(
    nInd = nFounders, # number of founder lines to be simulated
    nChr = 1, # number of chr
    segSites = nSnp, # number of snps
    inbred = TRUE, # if f.l. inbred
    species = "GENERIC" # all-purpose choice, historic bottlenecks
  )
# Define simulation parameters
SP <- SimParam$new(founderPop)
SP$setTrackPed(TRUE)
SP$setTrackRec(TRUE)
SP$addSnpChip(nSnpPerChr = nSnp)
# Input pedigree table and run population simulation
pop <- newPop(founderPop)
output <- pedigreeCross(pop,
                        id = pedigreeInitial[,2],
                        mother = pedigreeInitial[,4],
                        father = pedigreeInitial[,3])
# Get the output genotypes and pedigrees
geno <- pullSnpGeno(output)
pedigree <- data.frame(output@id,
                       output@father,
                       output@mother,
                       stringsAsFactors = FALSE)
masked_genotypes <- geno

## Add a 2% genotyping error. Can remove if errors do not need to be added. ONLY SIMULATE GES IN GENERATIONS G2-G4!!!!!

addErrors = function(values){
    newValues = (values + round(runif(length(values), min = 1, max = 2))) %% 3
    return(newValues)
}
G2onwards = masked_genotypes[firstG2:nrow(masked_genotypes),] ## NO GEs in G0G1
#markers = which(rbinom(length(G2onwards), 1, genErrorRate) == 1)
totalGEs = genErrorRate*length(G2onwards)
markers = sample(1:length(G2onwards), totalGEs , replace = FALSE) ## Randomly pick a proportion of the sites equal to error_rate for GEs
if(length(markers) >0){
    G2onwards[markers] = addErrors(G2onwards[markers])
}
masked_genotypes = rbind(masked_genotypes[1:(firstG2-1),], G2onwards)

# Calculates and stores the true recombination history. recHist contains the recombination event locations and recombinations contains the number of RE
recHist <- SP$recHist
countRecomb <- function(index) {
  ind = recHist[[index]]
  if(!is.null(ind)) {
    pat = ind[[1]][[2]]
    values_pat = nrow(pat) - 1 
    mat = ind[[1]][[1]]
    values_mat = nrow(mat) - 1 
    values_total = values_pat + values_mat
    return(c(index, values_pat, values_mat, values_total))
  }
  return(c(index, NA, NA, NA))
}
tmp = lapply((nFounders+1):length(recHist), countRecomb)
recombinations = do.call("rbind", tmp)
father <- c()
mother <- c()
gen <- c()
for (i in 1:nrow(recombinations)) {
  ind <- pedigreeInitial[recombinations[i, 1], 2]
  father <- append(father, pedigreeInitial[pedigreeInitial[, 2] %in% ind, 3])
  mother <- append(mother, pedigreeInitial[pedigreeInitial[, 2] %in% ind, 4])
  gen <- append(gen, pedigreeInitial[pedigreeInitial[, 2] %in% ind, 1])
}
recombinations_copy <- data.frame(gen, pedigreeInitial[recombinations[,1], 2], father, mother, recombinations[,2:4])

# Create a table showing the RE locations in every meiosis
RElocations <- c()
meiosis <- c()
indID <- c()
gen <- c()
for (i in (nFounders+1):nrow(pedigree)) {
  patRE <- recHist[[i]][[1]][[2]]
  matRE <- recHist[[i]][[1]][[1]]
  ifelse(!(nrow(patRE) == 1), RElocations <- append(RElocations, patRE[-1, 2] - 1), 0) # Paternal meiosis
  ifelse(!(nrow(matRE) == 1), RElocations <- append(RElocations, matRE[-1, 2] - 1), 0) # Maternal meiosis
  meiosis <- append(meiosis, rep('P', nrow(patRE) - 1))
  meiosis <- append(meiosis, rep('M', nrow(matRE) - 1))
  indID <- append(indID, rep(pedigree[i, 1], (nrow(patRE) - 1) + (nrow(matRE) - 1) ))
  gen <- append(gen, rep(pedigreeInitial[pedigreeInitial[, 2] %in% pedigree[i, 1], 1], (nrow(patRE) - 1) + (nrow(matRE) - 1) ))
}
RElocationsTable <- data.frame(indID, RElocations, meiosis)
father <- c()
mother <- c()
for (i in 1:nrow(RElocationsTable)) {
  ind <- RElocationsTable[i, 1]
  father <- append(father, pedigree[pedigree[, 1] %in% ind, 2])
  mother <- append(mother, pedigree[pedigree[, 1] %in% ind, 3])
}
RElocationsTable <- data.frame(gen, indID, father, mother, RElocations, meiosis)
# Write tables
for (col in 1:ncol(masked_genotypes)) {
  masked_genotypes[masked_genotypes[,col] == 0, col] = '1 1'
  masked_genotypes[masked_genotypes[,col] == 1, col] = '1 2'
  masked_genotypes[masked_genotypes[,col] == 2, col] = '2 2'
}
masked_genotypes = as.data.frame(masked_genotypes)
for (col in 1:ncol(masked_genotypes)) {
  masked_genotypes[,col] = stringr::str_split_fixed(masked_genotypes[,col], ' ', 2)
}
write.table(cbind(pedigreeInitial, data.frame(phenotype = 0, sex = 0), masked_genotypes),
            file = file.path(paste0(simID, ".ped")),
            sep = " ",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)
write.table(mapDF,
            file = file.path(paste0(simID, ".map")),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
write.table(RElocationsTable,
            file = file.path(paste0(simID, ".alpha.reco")),
            quote = FALSE,
            col.names = c('gen', 'id', 'ft', 'mt', 'RE', 'sex'),
            row.names = FALSE)
