#!/usr/bin/env Rscript

#=================================================================================================
#title: 6.recombination.frequency_function.R
#description: calculates recombination frequency from the recombination events detected and outputs .freq file based on one of the three methods: Poisson, cake and midpoint
#author: jmontero
#email: jose.a.montero-tena@ab.uni-giessen.de
#date: 2024-01-25
#version: 1.0.0
#notes: can be executed from 0.haplomagic.R
#	Output unit in cM/Kbp
#=================================================================================================

CalculateRF_bins = function(pop, chr, binSizeSNPs = 10, subset = 'RE') {
  recoDF = read.table(paste0(pop, "_", chr, ".f.reco"), sep = " ", header = TRUE, colClasses = c(rep("character", 11), rep("numeric", 4), rep("logical", 2)))
  mapDF = read.table(paste0(pop, "_", chr, ".map"), sep = "\t", header = FALSE, colClasses = c("character", "character", "numeric", "numeric"), col.names = c("chr", "snp", "gen", "pos"))
  pedDF = read.table(paste0(pop, "_", chr, ".pedigree"), sep = " ", header = TRUE, colClasses = "character")
  
  # Add the crossover midpoints
  recoDF$cMP = round((recoDF$cStart + recoDF$cEnd)/2, 0)
  
  # Subset the map keeping the adjacent markers of the bins
  mapDF = mapDF[c(seq(1, nrow(mapDF), binSizeSNPs), nrow(mapDF)),]
  
  # Add the interval distance to the map DF
  mapDF$interval.bp = c(mapDF$pos[2:length(mapDF$pos)] - mapDF$pos[1:length(mapDF$pos)-1], NA)
  
  # Create rec freq DF template
  freqDF = data.frame(pop = rep(pop, nrow(mapDF)),
                      chr = mapDF$chr[mapDF$chr == chr],
                      snp = mapDF$snp[mapDF$chr == chr],
                      aPos = 1:nrow(mapDF),
                      pPos = mapDF$pos[mapDF$chr == chr],
                      interval.bp = mapDF$interval.bp[mapDF$chr == chr])
  
  # Create matrix and fill them up to calculate recombination frequency in several steps
  if ( subset == "CO" ) { subset(recoDF, CO == TRUE) } else if ( subset == "GC" ) { recoDF = subset(recoDF, GC == TRUE) }
  
  #### INITIALIZATION
  ## Create count matrix {recombination event x snp interval}
  countMatrix = matrix(0,
                       nrow = nrow(recoDF),
                       ncol = length(freqDF$interval.bp[!is.na(freqDF$interval.bp)]))
  #### ITERATION
  for (i in seq(1, nrow(recoDF), 1)) {
    ## Add 1 to the bin where the crossover midpoint is
    countMatrix[i, which(freqDF$pPos == 
                           max(freqDF$pPos[freqDF$pPos < recoDF$cMP[i] & !is.na(freqDF$pPos)]))] = 1
  }
  ## Sum the crossovers by bins to calculate the number of crossovers per bin
  countVector = colSums(countMatrix)
  
  ## Total number of informative genotypes
  TotalNumberInformativeGenotypes = length(pedDF$gen[!(pedDF$gen %in% c("G0", "G1"))])*2 ## *2=Two meiosis per individual (P/M)
  
  ## Calculate recombination in cM/Kbp following the formulate: (100*crossovers/total_meiosis) / (basepairs/1000)
  RecFreq = ( 100 * countVector / TotalNumberInformativeGenotypes ) / ( freqDF$interval.bp[-length(freqDF$interval.bp)] * (1/1000) )
  
  ## Sum up to get individual x SNP information and return dataframe
  freqDF2 = data.frame(c(countVector, NA),
                       c(RecFreq, NA))
  ## Add info
  colnames(freqDF2) = c(paste0(subset, ".bins.noFilt"), paste0(subset, ".RecFreq.cM.Kbp.bins"))
  out = cbind(freqDF, freqDF2)
  out = out[!is.na(out$interval.bp), ]
  # Output the table with the RF for every REtype alongside the marker info
  write.table(out, file = paste0(pop, "_", chr, ".bins.freq"), sep = " ", quote = FALSE, col.names = TRUE, row.names = FALSE)
}

CalculateRF_Poisson = function(pop, chr, subset) {
  recoDF = read.table(paste0(pop, "_", chr, ".f.reco"), sep = " ", header = TRUE, colClasses = c(rep("character", 11), rep("numeric", 4), rep("logical", 2)))
  mapDF = read.table(paste0(pop, "_", chr, ".map"), sep = "\t", header = FALSE, colClasses = c("character", "character", "numeric", "numeric"), col.names = c("chr", "snp", "gen", "pos"))
  pedDF = read.table(paste0(pop, "_", chr, ".pedigree"), sep = " ", header = TRUE, colClasses = "character")
  
  # Add the interval distance to the map DF
  mapDF$interval.bp = c(mapDF$pos[2:length(mapDF$pos)] - mapDF$pos[1:length(mapDF$pos)-1], NA)
  
  # Create rec freq DF template
  freqDF = data.frame(pop = rep(pop, nrow(mapDF)),
                      chr = mapDF$chr[mapDF$chr == chr],
                      snp = mapDF$snp[mapDF$chr == chr],
                      aPos = 1:nrow(mapDF),
                      pPos = mapDF$pos[mapDF$chr == chr],
                      interval.bp = mapDF$interval.bp[mapDF$chr == chr])
  
  # Create matrix and fill them up to calculate recombination frequency in several steps
  if ( subset == "CO" ) { subset(recoDF, CO == TRUE) } else if ( subset == "GC" ) { recoDF = subset(recoDF, GC == TRUE) }
  
  #### INITIALIZATION
  ## Create matrices to fill with unambiguous intervals where COs were detected and ambiguous ones {recombination event x snp interval}
  breakpointIntervals = matrix(0,
                               nrow = nrow(recoDF),
                               ncol = length(freqDF$interval.bp[!is.na(freqDF$interval.bp)]))
  ambiguousIntervals = matrix(0,
                              nrow = nrow(recoDF),
                              ncol = length(freqDF$interval.bp[!is.na(freqDF$interval.bp)]))
  #### ITERATION
  for (i in seq(1, nrow(recoDF), 1)) {
    REinterval = recoDF$aEnd[i] - recoDF$aStart[i]
    ## Count the event as breakpoint interval or ambiguous interval
    if (REinterval == 1) {
      breakpointIntervals[i, seq(recoDF$aStart[i], recoDF$aEnd[i] - 1, 1)] = 1
    } else {
      ambiguousIntervals[i, seq(recoDF$aStart[i], recoDF$aEnd[i] - 1, 1)] = 1
    }
  }
  freqDF2 = data.frame(breakpointIntervals = c(colSums(breakpointIntervals), 0),
                       ambiguousIntervals = c(colSums(ambiguousIntervals), 0))
  totalNumberInformativeGenotypes = length(pedDF$gen[!(pedDF$gen %in% c("G0", "G1"))])*2 ## *2=Two meiosis per individual (P/M)
  freqDF2$totalNumberUnambiguousIntervals = totalNumberInformativeGenotypes - freqDF2$ambiguousIntervals
  freqDF2$negativeIntervals = freqDF2$totalNumberUnambiguousIntervals - freqDF2$breakpointIntervals
  freqDF2$averageNumberBreakpointsPerInterval = freqDF2$breakpointIntervals / freqDF2$negativeIntervals
  freqDF2$averageNumberCOsPerInterval = freqDF2$averageNumberBreakpointsPerInterval * freqDF2$totalNumberUnambiguousIntervals
  #freqDF2$cM = freqDF2$averageNumberCOsPerInterval * 100 / totalNumberInformativeGenotypes
  ## For cM calculation, divide by UNAMBIGUOUS GAMETES (NOT TOTAL GAMETES, BUT THOSE THAT CAN BE ASSIGNED TO A PARENT)
  freqDF2$cM = freqDF2$averageNumberCOsPerInterval * 100 / freqDF2$totalNumberUnambiguousIntervals
  freqDF2$RF = freqDF2$cM * 1000 / freqDF$interval.bp 
  keep = freqDF2[, c('averageNumberCOsPerInterval', 'RF')]
  colnames(keep) = c(paste0(subset, ".Poisson.noFilt"), paste0(subset, ".RecFreq.cM.Kbp.Poisson"))
  out = cbind(freqDF, keep)
  out = out[!is.na(out$interval.bp), ]
  # Output the table with the RF for every REtype alongside the marker info
  write.table(out, file = paste0(pop, "_", chr, ".Poisson.freq"), sep = " ", quote = FALSE, col.names = TRUE, row.names = FALSE)
}

CalculateRF_cake = function(pop, chr, subset) {
  recoDF = read.table(paste0(pop, "_", chr, ".f.reco"), sep = " ", header = TRUE, colClasses = c(rep("character", 11), rep("numeric", 4), rep("logical", 2)))
  mapDF = read.table(paste0(pop, "_", chr, ".map"), sep = "\t", header = FALSE, colClasses = c("character", "character", "numeric", "numeric"), col.names = c("chr", "snp", "gen", "pos"))
  pedDF = read.table(paste0(pop, "_", chr, ".pedigree"), sep = " ", header = TRUE, colClasses = "character")

  # Add the interval distance to the map DF
  mapDF$interval.bp = c(mapDF$pos[2:length(mapDF$pos)] - mapDF$pos[1:length(mapDF$pos)-1], NA)

    # Create rec freq DF template
  freqDF = data.frame(pop = rep(pop, nrow(mapDF)),
                      chr = mapDF$chr[mapDF$chr == chr],
                      snp = mapDF$snp[mapDF$chr == chr],
                      aPos = 1:nrow(mapDF),
                      pPos = mapDF$pos[mapDF$chr == chr],
                      interval.bp = mapDF$interval.bp[mapDF$chr == chr])

  # Create matrix and fill them up to calculate recombination frequency in several steps
  if ( subset == "CO" ) { subset(recoDF, CO == TRUE) } else if ( subset == "GC" ) { recoDF = subset(recoDF, GC == TRUE) }

  #### INITIALIZATION
  ## Create count matrix {recombination event x snp interval}
  countMatrix = matrix(0,
                       nrow = nrow(recoDF),
                       ncol = length(freqDF$interval.bp[!is.na(freqDF$interval.bp)]))
  ## Copy for factor matrix
  factorMatrix = countMatrix
  #### ITERATION
  for (i in seq(1, nrow(recoDF), 1)) {
    REinterval = recoDF$aEnd[i] - recoDF$aStart[i]
    ## Fill up the count matrix with "1/number of SNPs in RE interval" ()
    countMatrix[i, seq(recoDF$aStart[i], recoDF$aEnd[i] - 1, 1)] =
      1/REinterval
    ## Fill up the factor matrix with interval size and number of genotypes
    TotalNumberInformativeGenotypes = length(pedDF$gen[!(pedDF$gen %in% c("G0", "G1"))])*2 ## *2=Two meiosis per individual (P/M)
    REintervalDistance = recoDF$cEnd[i] - recoDF$cStart[i]
    factorMatrix[i, seq(recoDF$aStart[i], recoDF$aEnd[i] - 1, 1)] =
      100 / (REintervalDistance / 1e3 * TotalNumberInformativeGenotypes) ## 100/ = Convert from Morgan to centiMorgan
  }
  #### CONVERGENCE
  ## no.recombinants * factor (1/total.ids x distance) = RF0    
  RecFreqMatrix = countMatrix * factorMatrix
  ## Sum up to get individual x SNP information and return dataframe
  freqDF2 = data.frame(c(colSums(countMatrix), NA),
                       c(colSums(RecFreqMatrix), NA))
  colnames(freqDF2) = c(paste0(subset, ".Cake.noFilt"), paste0(subset, ".RecFreq.cM.Kbp.Cake"))
  out = cbind(freqDF, freqDF2)
  out = out[!is.na(out$interval.bp), ]
  # Output the table with the RF for every REtype alongside the marker info
  write.table(out, file = paste0(pop, "_", chr, ".Cake.freq"), sep = " ", quote = FALSE, col.names = TRUE, row.names = FALSE)
}

CalculateRF_mp = function(pop, chr, subset) {
  recoDF = read.table(paste0(pop, "_", chr, ".f.reco"), sep = " ", header = TRUE, colClasses = c(rep("character", 11), rep("numeric", 4), rep("logical", 2)))
  mapDF = read.table(paste0(pop, "_", chr, ".map"), sep = "\t", header = FALSE, colClasses = c("character", "character", "numeric", "numeric"), col.names = c("chr", "snp", "gen", "pos"))
  pedDF = read.table(paste0(pop, "_", chr, ".pedigree"), sep = " ", header = TRUE, colClasses = "character")

  # Add the interval distance to the map DF
  mapDF$interval.bp = c(mapDF$pos[2:length(mapDF$pos)] - mapDF$pos[1:length(mapDF$pos)-1], NA)

    # Create rec freq DF template
  freqDF = data.frame(pop = rep(pop, nrow(mapDF)),
                      chr = mapDF$chr[mapDF$chr == chr],
                      snp = mapDF$snp[mapDF$chr == chr],
                      aPos = 1:nrow(mapDF),
                      pPos = mapDF$pos[mapDF$chr == chr],
                      interval.bp = mapDF$interval.bp[mapDF$chr == chr])

  # Create matrix and fill them up to calculate recombination frequency in several steps
  if ( subset == "CO" ) { subset(recoDF, CO == TRUE) } else if ( subset == "GC" ) { recoDF = subset(recoDF, GC == TRUE) }


  #### INITIALIZATION
  ## Create count matrix {recombination event x snp interval}
  countMatrix = matrix(0,
                       nrow = nrow(recoDF),
                       ncol = length(freqDF$interval.bp[!is.na(freqDF$interval.bp)]))
  ## Copy for factor matrix
  factorMatrix = countMatrix
  #### ITERATION
  for (i in seq(1, nrow(recoDF), 1)) {
    REinterval = 1
    ## Fill up the count matrix with "1/number of SNPs in RE interval" ()
    countMatrix[i, round(mean(c(recoDF$aStart[i], recoDF$aEnd[i]-1)), 0)] =
      1/REinterval
    ## Fill up the factor matrix with interval size and number of genotypes
    TotalNumberInformativeGenotypes = length(pedDF$gen[!(pedDF$gen %in% c("G0", "G1"))])*2 ## *2=Two meiosis per individual (P/M)
    REintervalDistance = mapDF$pos[round(mean(c(recoDF$aStart[i], recoDF$aEnd[i]-1)))+1] - mapDF$pos[round(mean(c(recoDF$aStart[i], recoDF$aEnd[i]-1)))]
    factorMatrix[i, round(mean(c(recoDF$aStart[i], recoDF$aEnd[i]-1)), 0)] =
      100/( REintervalDistance / 1e3 * TotalNumberInformativeGenotypes ) ## 100/ = Convert from Morgan to centiMorgan
  }
  #### CONVERGENCE
  ## no.recombinants * factor (1/total.ids x distance) = RF
  RecFreqMatrix = countMatrix * factorMatrix
  ## Sum up to get individual x SNP information and return dataframe
  freqDF2 = data.frame(c(colSums(countMatrix), NA),
                       c(colSums(RecFreqMatrix), NA))
  colnames(freqDF2) = c(paste0(subset, ".MP.noFilt"), paste0(subset, ".RecFreq.cM.Kbp.MP"))
  out = cbind(freqDF, freqDF2)
  out = out[!is.na(out$interval.bp), ]
  # Output the table with the RF for every REtype alongside the marker info
  write.table(out, file = paste0(pop, "_", chr, ".MP.freq"), sep = " ", quote = FALSE, col.names = TRUE, row.names = FALSE)
}
