#!/usr/bin/env Rscript

#=================================================================================================
#title: 5.statistics_function.R
#description: calculates statistics from the recombination events detected and outputs .stats file
#author: jmontero
#email: jose.a.montero-tena@ab.uni-giessen.de
#date: 2023-02-21
#version: 1.0.0
#notes: can be executed from 0.haplomagic.R
#=================================================================================================

GetHaplomagicStatistics = function(pop, chr) {
  # Import relevant files
  pedDF = read.table(paste0(pop, "_", chr, ".pedigree"), sep = " ", header = TRUE, colClasses = "character")
  originMatrix.noFilt.noImp = as.matrix(read.table(paste0(pop, "_", chr, ".nf.ni.origin"), sep = " ", header = FALSE, colClasses = "character", row.names = 1))
  recoDF = read.table(paste0(pop, "_", chr, ".f.reco"), sep = " ", header = TRUE, colClasses = c(rep("character", 11), rep("numeric", 4), rep("logical", 2)))
  recoDF.noFilt = read.table(paste0(pop, "_", chr, ".nf.reco"), sep = " ", header = TRUE, colClasses = c(rep("character", 11), rep("numeric", 4), rep("logical", 2)))
  nSnp = ncol(originMatrix.noFilt.noImp) # Shortcut to number of SNPs
  
  statsDF = data.frame()
  for (i in seq(1, nrow(originMatrix.noFilt.noImp), 1)) {
    id = strsplit(rownames(originMatrix.noFilt.noImp)[i], "_")[[1]][1]
    gen = pedDF[pedDF$id == id, "gen"]
    sire = pedDF[pedDF$id == id, "sire"]
    dam = pedDF[pedDF$id == id, "dam"]
    mei = strsplit(rownames(originMatrix.noFilt.noImp)[i], "_")[[1]][2]
    RE.noFilt = nrow(recoDF.noFilt[recoDF.noFilt$id == id & recoDF.noFilt$mei == mei,])
    CO.noFilt = nrow(recoDF.noFilt[recoDF.noFilt$id == id & recoDF.noFilt$mei == mei & recoDF.noFilt$CO == TRUE,])
    GC.noFilt = nrow(recoDF.noFilt[recoDF.noFilt$id == id & recoDF.noFilt$mei == mei & recoDF.noFilt$GC == TRUE,])
    RE.filt = nrow(recoDF[recoDF$id == id & recoDF$mei == mei,])
    CO.filt = nrow(recoDF[recoDF$id == id & recoDF$mei == mei & recoDF$CO == TRUE,])
    GC.filt = nrow(recoDF[recoDF$id == id & recoDF$mei == mei & recoDF$GC == TRUE,])
    statsDFTMP = data.frame(
      pop = pop,
      chr = chr,
      gen = gen,
      id = id,
      sire = sire,
      dam = dam,
      mei = mei,
      RE.noFilt = RE.noFilt,
      CO.noFilt = CO.noFilt,
      GC.noFilt = GC.noFilt,
      RE.filt = RE.filt,
      CO.filt = CO.filt,
      GC.filt = GC.filt,
      RE.fRate = ifelse(RE.noFilt == 0,
                        NA,
                        round(1 - RE.filt/RE.noFilt, 2)),
      CO.fRate = ifelse(CO.noFilt == 0,
                        NA,
                        round(1 - CO.filt/CO.noFilt, 2)),
      GC.fRate = ifelse(GC.noFilt == 0,
                        NA,
                        round(1 - GC.filt/GC.noFilt, 2)),
      nonInfRate = round(length(which(originMatrix.noFilt.noImp[i, ] == "*"))/nSnp, 2),
      MDrate = round(length(which(originMatrix.noFilt.noImp[i, ] == "?"))/nSnp, 2),
      MErate = round(length(which(originMatrix.noFilt.noImp[i, ] == "!"))/nSnp, 2)
    )
    statsDF = rbind(statsDF, statsDFTMP)
  }
  write.table(statsDF, file = paste0(pop, "_", chr, ".stats"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep =  " ")
}

