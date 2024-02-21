#!/usr/bin/env Rscript

#================================================================================
#title: 4.detect.recombination_function.R
#description: detects recombination events from .origin file & outputs .reco file
#author: jmontero
#email: jose.a.montero-tena@ab.uni-giessen.de
#date: 2023-02-21
#version: 1.0.0
#notes: can be executed from 0.haplomagic.R
#================================================================================

DetectRecombination = function(pop, chr, thr) {
  # Import relevant files
  mapDF = read.table(paste0(pop, "_", chr, ".map"), sep = "\t", header = FALSE, colClasses = c("character", "character", "numeric", "numeric"), col.names = c("chr", "id", "gen", "phy"))
  pedDF = read.table(paste0(pop, "_", chr, ".pedigree"), sep = " ", header = TRUE, colClasses = "character")
  originMatrix = as.matrix(read.table(paste0(pop, "_", chr, ".f.i.origin"), sep = " ", header = FALSE, colClasses = "character", row.names = 1))
  originMatrix.noFilt = as.matrix(read.table(paste0(pop, "_", chr, ".nf.i.origin"), sep = " ", header = FALSE, colClasses = "character", row.names = 1))
  haploMatrix = as.matrix(read.table(paste0(pop, "_", chr, ".haplo"), sep = " ", header = FALSE, colClasses = "character", row.names = 1))
  nSnp = nrow(mapDF) # Shortcut to number of SNPs
  
  # Create list with filtered and nonfiltered data
  origin = list(originMatrix, originMatrix.noFilt)
  outName = c(paste0(pop, "_", chr, ".f.reco"), paste0(pop, "_", chr, ".nf.reco"))
  for (j in seq(1, length(origin), 1)) {
    # Create main dataframe from which the file will be created in the end
    recoDF = data.frame()
    # Iterate over the rows of the haplotypes of >G2 individuals
    for (i in 1:nrow(origin[[j]])) {
      # First get a vector with the start and end coordinates of the haplotype blocks
      getCoords = function (block) {
        return(ifelse(c(block[1] == 1, block[length(block)] == nSnp),
                      c(NA, NA),
                      block[c(1, length(block))]))
      }
      haplotypeIndexP = which(origin[[j]][i, ] == "P")
      haplotypeIndexP = split(haplotypeIndexP, cumsum(c(1, diff(haplotypeIndexP) != 1)))
      haplotypeIndexP = sort(unlist(lapply(haplotypeIndexP, getCoords)), decreasing = FALSE)
      haplotypeIndexM = which(origin[[j]][i, ] == "M")
      haplotypeIndexM = split(haplotypeIndexM, cumsum(c(1, diff(haplotypeIndexM) != 1)))
      haplotypeIndexM = sort(unlist(lapply(haplotypeIndexM, getCoords)), decreasing = FALSE)
      recCoords = sort(c(haplotypeIndexP[!is.na(haplotypeIndexP)], haplotypeIndexM[!is.na(haplotypeIndexM)]), decreasing = FALSE)
      if (length(recCoords) > 1) {
        # Save the array and chromosome positions
        recArrayStart = recCoords[seq(1, length(recCoords), 2)]
        recArrayEnd = recCoords[seq(2, length(recCoords), 2)]
        recChrStart = mapDF[recArrayStart, "phy"]
        recChrEnd = mapDF[recArrayEnd, "phy"]
        # Classify RE as CO/GC according to the provided threshold
        if (length(recCoords) == 2) {
          isHaploblockGC = c(FALSE, FALSE) # If there is only one RE, therefore two haploblocks, they will be classified as CO
        } else {
          isHaploblockGC = ifelse(recChrStart[seq(2, length(recChrStart), 1)] - recChrStart[seq(1, length(recChrStart)-1, 1)] < thr,
                                  TRUE, # If the haploblock length or distance between one RE start and the last RE end is smaller than the threshold, then TRUE GC
                                  FALSE)
          isHaploblockGC = c(FALSE, isHaploblockGC, FALSE) # The first and last haploblocks, just flanked by one RE, will be considered as CO
        }
        GC = ifelse(isHaploblockGC[seq(2, length(isHaploblockGC), 1)] | isHaploblockGC[seq(1, length(isHaploblockGC)-1, 1)],
                    TRUE, # If any of the flanking haploblocks is GC, the RE will be classified as GC
                    FALSE)
        CO = ifelse(isHaploblockGC[seq(2, length(isHaploblockGC), 1)] & isHaploblockGC[seq(1, length(isHaploblockGC)-1, 1)],
                    FALSE, 
                    TRUE) # If any of the flanking haploblocks is CO, the RE will be counted as CO. This means that RE can be classified as CO and GC if they are surrounded by either
        # Save information about where the RE happened
        id = strsplit(rownames(origin[[j]])[i], "_")[[1]][1]
        gen = pedDF[pedDF$id == id, "gen"]
        sire = pedDF[pedDF$id == id, "sire"]
        dam = pedDF[pedDF$id == id, "dam"]
        meiosis = strsplit(rownames(origin[[j]])[i], "_")[[1]][2]
        origin1 = origin[[j]][i, recArrayStart]
        origin2 = origin[[j]][i, recArrayEnd]
        haplo1 = haploMatrix[paste0(id, "_", meiosis), recArrayStart]
        haplo2 = haploMatrix[paste0(id, "_", meiosis), recArrayEnd]
        # Create dataframe and bind it to the main dataframe
        chrDF = data.frame(pop = rep(pop, length(recArrayStart)),
                           chr = rep(chr, length(recArrayStart)),
                           gen = rep(gen, length(recArrayStart)),
                           id = rep(id, length(recArrayStart)),
                           sire = rep(sire, length(recArrayStart)),
                           dam = rep(dam, length(recArrayStart)),
                           mei = rep(meiosis, length(recArrayStart)),
                           ori1 = origin1,
                           ori2 = origin2,
                           haplo1 = haplo1,
                           haplo2 = haplo2,
                           aStart = recArrayStart,
                           aEnd = recArrayEnd,
                           cStart = recChrStart,
                           cEnd = recChrEnd,
                           CO = CO,
                           GC = GC)
        recoDF = rbind(recoDF, chrDF)
      }
    }
    # Write table
    write.table(recoDF, paste0(outName[j]), sep = " ", row.names = FALSE, quote = FALSE, col.names = TRUE)
  }
}
