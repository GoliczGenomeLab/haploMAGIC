#!/usr/bin/env Rscript

#=======================================================================================================================================
#title: get.performance.parameters.R
#description: calculates precision, recall, median recombination interval size and proportion of RE removed in the simulated populations
#author: jmontero
#date: 2023-03-20
#version: 1.0.0
#usage: Rscript 2.phase.rest.R "pop1 pop2 pop3" "chr1 chr2 chr3" 2 None noCor
#notes: requires the input the arguments min, imp and cor
#	assumes that the simulated population name is simpop_GErate_population_chr
#=======================================================================================================================================

args = commandArgs(trailingOnly=TRUE)

# Define global variables
pops=unlist(strsplit(args[1], split = " "))
chrs=as.numeric(unlist(strsplit(args[2], split = " ")))
min=2
imp=None
cor=noCor

for (pop in pops) {
  for (chr in chrs) {
    # Upload the tables with the RecoPed and the simulation data
    if (min == 1) {
      reco = read.table(paste0(pop, "_", chr, ".nf.reco"), colClasses = c(rep("character", 11), rep("numeric", 4), rep("logical", 2)), header = TRUE)
    } else if (min > 1) {
      reco = read.table(paste0(pop, "_", chr, ".f.reco"), colClasses = c(rep("character", 11), rep("numeric", 4), rep("logical", 2)), header = TRUE)
    }
    nfreco = read.table(paste0(pop, "_", chr, ".nf.reco"), colClasses = c(rep("character", 11), rep("numeric", 4), rep("logical", 2)), header = TRUE)
    reco$cLength = reco$cEnd - reco$cStart
    reco$aLength = reco$aEnd - reco$aStart
    statsDF = read.table(paste0(pop, "_", chr, ".stats"), header = TRUE, colClasses = c(rep("character", 5), rep("numeric", 12)))
    statsDF$infRate = statsDF$nonInfRate + statsDF$MDrate + statsDF$MErate
    alpha = read.table(paste0(pop, "_", chr, ".alpha.reco"), colClasses = c(rep("character", 4), rep("numeric", 1), rep("character", 1)), header = TRUE)
    alpha = alpha[alpha$gen != "G1", ]
    
    # Create empty vectors precision and recall to fill up with every individual's data
    
    precision = c()
    recall = c()
    
    gen = c()
    ids = c()
    sires = c()
    dams = c()
    mei = c()
    Dglobal = c()
    Sglobal = c()
    SDglobal = c()
    
    # Iterate over the individuals and compare the simulated and detected events
    individualList = unique(statsDF$id)
    for (id in individualList) {
      # Create the subsets with the rows of the individuals
      for (meiosis in c("P", "M")) {
        
        gen = append(gen, statsDF[statsDF$id == id, "gen"][1])
        sires = append(sires, statsDF[statsDF$id == id, "sire"][1])
        dams = append(dams, statsDF[statsDF$id == id, "dam"][1])

        ids = append(ids, id)
        mei = append(mei, meiosis)
        
        recoSub = reco[reco$id == id & reco$mei == meiosis,]
        alphaSub = alpha[alpha$id == id & alpha$sex == meiosis,]
        # Count D, number of detected events, and S, simulated
        D = nrow(recoSub)
        S = nrow(alphaSub)
        SD = 0
        Dglobal = append(Dglobal, D)
        Sglobal = append(Sglobal, S)
        SD = 0 # NUMBER OF SIMULATED EVENTS DETECTED, THAT WILL BE FILLED BY ITERATIONS IF THE DETECTED VALUES ARE WITHIN THE SIMULATED ONES
        # Iterate over simulated events
        if (D > 0 & S > 0) {
          discard = c()
          for (sim in alphaSub$RE) {
            index = seq(1, D)
            index = index[!index %in% discard] # Eliminate the index of this matching row for avoiding rematching and boosting SD
            for (i in index) {
              if (sim >= recoSub$aStart[i] & sim < recoSub$aEnd[i]) {
                SD = SD + 1
                discard = append(discard, i)
                break
              }
            }
          }
        }
        SDglobal = append(SDglobal, SD)
      }
    }
    
    if (sum(Sglobal) == 0) {
      recall = NA
    } else {
      recall = sum(SDglobal) / sum(Sglobal)
    }
    
    if (sum(Dglobal) == 0) {
      precision = NA
    } else {
      precision = sum(SDglobal) / sum(Dglobal)
    }
    F1score = 2*precision*recall/(precision+recall)
    GErate = strsplit(pop, '_')[[1]][2]
    medianREintervalSize = median(reco$aLength)
    filteredREs = 1-nrow(reco)/nrow(nfreco)
    cat(paste0("haploMAGIC ", GErate, " ", pop, " ", chr, " ", min, " ", imp, " ", cor, " ", precision, " ", recall, " ", F1score, " ", medianREintervalSize, " ",  filteredREs, "\n"))
  }
}
