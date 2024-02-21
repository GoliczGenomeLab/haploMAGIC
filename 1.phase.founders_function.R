#!/usr/bin/env Rscript

#====================================================================================================
#title: 1.phase.founders_function.R
#description: phases and haplotypes G0 and G1 individuals & creates files for future haploMAGIC steps
#author: jmontero
#email: jose.a.montero-tena@ab.uni-giessen.de
#date: 2023-02-21
#version: 1.0.0
#notes: can be executed from 0.haplomagic.R
#	.haplo file contains letters representing G0 founders
#====================================================================================================

PhaseFounders = function(pop, chr) {
    # Import PED file
    pedFileName = paste0(pop, "_", chr, ".ped")
    pedFile = read.table(pedFileName, sep = " ", header = FALSE, colClasses = "character")
    nSnp = (ncol(pedFile) - 6)/2 # Shortcut to number of SNPs
    
    # Create a pedigree data frame
    pedDF = pedFile[, 1:4]
    colnames(pedDF) = c("gen", "id", "sire", "dam")
    
    # Create a genotype matrix
    genMatrix = as.matrix(pedFile[, 7:ncol(pedFile)])
    rownames(genMatrix) = pedDF$id
    
    # Create a phase matrix
    phaseMatrix = matrix(nrow = 0, ncol = nSnp)
    founderLines = pedDF[pedDF$gen == "G0", "id"]
    founderRowNames = c()
    for (founder in founderLines) {
      patPhase = genMatrix[founder, seq(1, ncol(genMatrix), 2)]
      matPhase = genMatrix[founder, seq(2, ncol(genMatrix), 2)]
      if (all(patPhase == matPhase)) { # Check if founders are homozygous
        phaseMatrix = rbind(phaseMatrix, 
                            patPhase, 
                            matPhase)
        founderRowNames = append(founderRowNames, paste0(founder, c("_P", "_M")))
      } else {
        cat(paste0("Error: Founder line ", founder, " is not homozygous"))
      }
    }
    rownames(phaseMatrix) = founderRowNames
    
    # Add G1 to the phase matrix
    phaseMatrixTMP = matrix(nrow = 0, ncol = nSnp)
    G1ids = pedDF[pedDF$gen == "G1", "id"]
    G1RowNames = c()
    for (id in G1ids) {
      sire = pedDF[pedDF$id == id, "sire"]
      dam = pedDF[pedDF$id == id, "dam"]
      patPhase = phaseMatrix[paste0(sire, "_P"), ]
      matPhase = phaseMatrix[paste0(dam, "_P"), ]
      phaseMatrixTMP = rbind(phaseMatrixTMP, 
                             patPhase, 
                             matPhase)
      G1RowNames = append(G1RowNames, paste0(id, c("_P", "_M")))
    }
    rownames(phaseMatrixTMP) = G1RowNames
    phaseMatrix = rbind(phaseMatrix, phaseMatrixTMP)
    
    # Create haploblock matrix
    labels = c(toupper(letters), letters, 0:1000)
    haploMatrix = matrix(nrow = 0, ncol = nSnp)
    founderLines = pedDF[pedDF$gen == "G0", "id"]
    for (f in seq(1, length(founderLines), 1)) {
      labelSeq = rep(labels[f], nSnp)
      haploMatrix = rbind(haploMatrix, 
                          labelSeq,
                          labelSeq)
    }
    rownames(haploMatrix) = founderRowNames
    
    # Add G1 to haploblock matrix
    haploMatrixTMP = matrix(nrow = 0, ncol = nSnp)
    for (id in G1ids) {
      sire = pedDF[pedDF$id == id, "sire"]
      dam = pedDF[pedDF$id == id, "dam"]
      patHaplo = haploMatrix[paste0(sire, "_P"), ]
      matHaplo = haploMatrix[paste0(dam, "_P"), ]
      haploMatrixTMP = rbind(haploMatrixTMP, 
                             patHaplo, 
                             matHaplo)
    }
    rownames(haploMatrixTMP) = G1RowNames
    haploMatrix = rbind(haploMatrix, haploMatrixTMP)
    
    write.table(pedDF, paste0(pop, "_", chr, ".pedigree"), sep = " ", row.names = FALSE, quote = FALSE, col.names = TRUE)
    write.table(genMatrix, paste0(pop, "_", chr, ".genotype"), sep = " ", row.names = TRUE, quote = FALSE, col.names = FALSE)
    write.table(phaseMatrix, paste0(pop, "_", chr, ".phase"), sep = " ", row.names = TRUE, quote = FALSE, col.names = FALSE)
    write.table(haploMatrix, paste0(pop, "_", chr, ".haplo"), sep = " ", row.names = TRUE, quote = FALSE, col.names = FALSE)
}
