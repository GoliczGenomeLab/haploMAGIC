#!/usr/bin/env Rscript

#========================================================================================================
#title: 3.haploblock_function.R
#description: assigns founder haplotypes to every allele in the population
#author: jmontero
#email: jose.a.montero-tena@ab.uni-giessen.de
#date: 2023-02-21
#version: 1.0.0
#notes: can be executed from 0.haplomagic.R
#	.haplo contains letters representing G0 founders
#       no gap imputation (* in the parent => * in the child). Might result in big gaps in last generation
#=========================================================================================================

AssignFounderHaploblocks = function(pop, chr) {
  # Import relevant files
  pedDF = read.table(paste0(pop, "_", chr, ".pedigree"), sep = " ", header = TRUE, colClasses = "character")
  haploMatrix = as.matrix(read.table(paste0(pop, "_", chr, ".haplo"), sep = " ", header = FALSE, colClasses = "character", row.names = 1))
  originMatrix = as.matrix(read.table(paste0(pop, "_", chr, ".f.i.origin"), sep = " ", header = FALSE, colClasses = "character", row.names = 1))
  nSnp = ncol(originMatrix) # Shortcut to number of SNPs
  
  # Iterate over generations to add rows to the haplomatrix from every individual
  for (gen in unique(pedDF$gen[!pedDF$gen %in% c('G0', 'G1')])) {
    ids = pedDF[pedDF$gen == gen, "id"]
    for (id in ids) {
      for (meiosis in c("P", "M")) {
        parent = ifelse(meiosis == "P",
                        pedDF[pedDF$id == id, "sire"],
                        pedDF[pedDF$id == id, "dam"])
        haploMatrixTMP = ifelse(originMatrix[paste0(id, "_", meiosis), ] == "*",
                                "*",
                                ifelse(originMatrix[paste0(id, "_", meiosis), ] == "P",
                                       haploMatrix[paste0(parent, "_P"), ],
                                       haploMatrix[paste0(parent, "_M"), ]))
        haploMatrixTMP = matrix(haploMatrixTMP, nrow = 1, ncol = nSnp)
        rownames(haploMatrixTMP) = paste0(id, "_", meiosis)
        haploMatrix = rbind(haploMatrix, haploMatrixTMP)
      }
    }
  }
  write.table(haploMatrix, paste0(pop, "_", chr, ".haplo"), sep = " ", row.names = TRUE, quote = FALSE, col.names = FALSE)
}
