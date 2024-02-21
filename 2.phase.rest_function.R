#!/usr/bin/env Rscript

#==========================================================================================
#title: 2.phase.rest_function.R
#description: phases and haplotypes individuals after G1 using phase and haplotypes from G1
#author: jmontero
#email: jose.a.montero-tena@ab.uni-giessen.de
#date: 2023-02-21
#version: 1.0.0
#notes: can be executed from 0.haplomagic.R
#	border filtering applied
#==========================================================================================

PhaseRest = function(pop, chr, min, imp, cor) {
  # Import relevant files
  pedDF = read.table(paste0(pop, "_", chr, ".pedigree"), sep = " ", header = TRUE, colClasses = "character")
  genMatrix = as.matrix(read.table(paste0(pop, "_", chr, ".genotype"), sep = " ", header = FALSE, colClasses = "character", row.names = 1))
  phaseMatrix = as.matrix(read.table(paste0(pop, "_", chr, ".phase"), sep = " ", header = FALSE, colClasses = "character", row.names = 1))
  haploMatrix = as.matrix(read.table(paste0(pop, "_", chr, ".haplo"), sep = " ", header = FALSE, colClasses = "character", row.names = 1))
  nSnp = ncol(phaseMatrix) # Shortcut to number of SNPs
  
  # Start by creating the haplotype origin matrix
  originMatrix = matrix(nrow = 0, ncol = nSnp) # Final matrix with imputed filtered haploblocks
  originMatrix.noFilt.noImp = originMatrix # Non-imputed non-filtered haploblocks for statistics
  originMatrix.noFilt = originMatrix # Non-filtered haploblocks for statistics
  for (gen in unique(pedDF$gen[!pedDF$gen %in% c('G0', 'G1')])) {
    ids = pedDF[pedDF$gen == gen, "id"]
    originMatrixTMP = matrix(nrow = 0, ncol = nSnp) ; phaseMatrixTMP = originMatrixTMP ; originMatrix.noFilt.noImpTMP = originMatrixTMP ; originMatrix.noFiltTMP = originMatrixTMP
    rowNames = c()
    for (id in ids) {
      # 1. Extract the id's genotype
      Igen = genMatrix[id, ]
      # 2. Extract the parents names
      sire = pedDF[pedDF$id == id, "sire"]
      dam = pedDF[pedDF$id == id, "dam"]
      # 3. Use the parents' names to extract their paternal and maternal phases
      FPphase = phaseMatrix[paste0(sire, "_P"), ]
      FMphase = phaseMatrix[paste0(sire, "_M"), ]
      MPphase = phaseMatrix[paste0(dam, "_P"), ]
      MMphase = phaseMatrix[paste0(dam, "_M"), ]
      # 4. PHASING
      # 4.1 Get the unphased child's alleles
      allele1 = Igen[seq(1, length(Igen), 2)] ; allele2 = Igen[seq(2, length(Igen), 2)] 
      # 4.2 Define phasing scenarios
      noMissingData = allele1 != 0 & allele2 != 0 & FPphase != 0 & FMphase != 0 & MPphase != 0 & MMphase != 0
noMendelianError = ((allele1 == FPphase | allele1 == FMphase) & (allele2 == MPphase | allele2 == MMphase)) | ((allele2 == FPphase | allele2 == FMphase) & (allele1 == MPphase | allele1 == MMphase)) # FIXED 23.03.14
      noTripleHet = !(allele1 != allele2 & FPphase != FMphase & MPphase != MMphase) # Triple heterozygous: Non-phaseable scenario in which all trio individuals are het
      homozygousID = allele1 == allele2
      homozygousFather = FPphase == FMphase
      homozygousMother = MPphase == MMphase
      # 4.3 Iterate over the ind phase to fill up the phase following the phasing scenarios
      IPphase = ifelse(noMissingData & noMendelianError, 
                       ifelse(homozygousID, 
                              allele1,
                              ifelse(homozygousFather, # If the mother is hom the father phase can be resolved
                                     FPphase,
                                     ifelse(homozygousMother, 
                                            ifelse(allele1 == MPphase, 
                                                   allele2,
                                                   allele1),
                                            "0"))),
                       "0")
      IMphase = ifelse(noMissingData & noMendelianError, 
                       ifelse(homozygousID, 
                                allele1,
                              ifelse(homozygousMother, # If the mother is hom the father phase can be resolved
                                     MPphase,
                                     ifelse(homozygousFather, 
                                            ifelse(allele1 == FPphase, 
                                                   allele2,
                                                   allele1),
                                            "0"))),
                       "0")
      # 5. HAPLOTYPING
      UnknownChildPaternalAllele = IPphase == "0" ; UnknownChildMaternalAllele = IMphase == 0
      IPorigin = ifelse(noMissingData & noMendelianError, 
                        ifelse(homozygousFather | UnknownChildPaternalAllele,
                               "*",
                               ifelse(IPphase == FPphase,
                                      "P",
                                      "M")),
                        ifelse(!noMissingData,
                               "?",
                               "!"))
      IMorigin = ifelse(noMissingData & noMendelianError, 
                        ifelse(homozygousMother | UnknownChildMaternalAllele,
                               "*",
                               ifelse(IMphase == MPphase,
                                      "P",
                                      "M")),
                        ifelse(!noMissingData,
                               "?",
                               "!"))
      # 6. HAPLOTYPE IMPUTATION
      Iorigin = list(IPorigin, IMorigin) ; Iphase = list(IPphase, IMphase) ; Fphase = list(FPphase, MPphase) ; Mphase = list(FMphase, MMphase)
      # For homozygous parents, avoid blank origin files by adding a random P
      for (i in seq(1, 2, 1)) {
        if (length(which(Iorigin[[i]] %in% c("P", "M"))) == 0) {
          Iorigin[[i]][1] = "P"
        }
      }
      for (i in seq(1, 2, 1)) {
        # Add the non-filtered non-imputed sequence to respective matrix
        originMatrix.noFilt.noImpTMP = rbind(originMatrix.noFilt.noImpTMP, Iorigin[[i]])
        # 6.1 1ST haplotype imputation round: impute non assigned alleles if they do not lie within RE gaps. ALSO save the indexes of the informative alleles
        wildcardIndex = which(Iorigin[[i]] == "*" | Iorigin[[i]] == "!" | Iorigin[[i]] == "?")
        informativeIndex = which(Iorigin[[i]] != "*" & Iorigin[[i]] != "!" & Iorigin[[i]] != "?") # This will be used to validate haploblocks based on the number of inf. alleles
        gaps = split(wildcardIndex, cumsum(c(1, diff(wildcardIndex) != 1)))
        if (length(gaps[[1]]) != 0 ) {
          for (gap in gaps) {
            Iorigin[[i]][seq(gap[1], gap[length(gap)], 1)] = ifelse(gap[1] == 1,
                                                                    Iorigin[[i]][gap[length(gap)] + 1],
                                                                    ifelse(gap[length(gap)] == nSnp,
                                                                           Iorigin[[i]][gap[1] - 1],
                                                                           ifelse(Iorigin[[i]][gap[1] - 1] == Iorigin[[i]][gap[length(gap)] + 1],
                                                                                  Iorigin[[i]][gap[length(gap)] + 1],
                                                                                  "*")))
          }
        } else {
          Iorigin[[i]] = Iorigin[[i]] # If there is no haplotype to impute
        }
        originMatrix.noFiltTMP = rbind(originMatrix.noFiltTMP, Iorigin[[i]])
        # 6.2 2ND haplotype imputation round: filter out haplotypes with no more than the minimum number of informative alleles (FILTERING)
        haplotypeIndexP = which(Iorigin[[i]] == "P") # This will be used to validate haploblocks based on the number of inf. alleles
        haplotypeIndexM = which(Iorigin[[i]] == "M")
        haploblocks = append(split(haplotypeIndexP, cumsum(c(1, diff(haplotypeIndexP) != 1))), split(haplotypeIndexM, cumsum(c(1, diff(haplotypeIndexM) != 1))))
        if (length(haploblocks[[1]]) != 0 ) {
          for (haploblock in haploblocks) {
            if (length(haploblock) != 0) {
              Iorigin[[i]][seq(haploblock[1], haploblock[length(haploblock)], 1)] = ifelse(length(informativeIndex[informativeIndex >= haploblock[1] & informativeIndex <= haploblock[length(haploblock)]]) < min,
                                                                                           "*",
                                                                                           Iorigin[[i]][seq(haploblock[1], haploblock[length(haploblock)], 1)])
            }
          }
        }

        # 6.3 3RD haplotype imputation round: Repetition of the 1st round to impute the filtered haplotypes
        wildcardIndex = which(Iorigin[[i]] == "*" | Iorigin[[i]] == "!" | Iorigin[[i]] == "?")
        informativeIndex = which(Iorigin[[i]] != "*" & Iorigin[[i]] != "!" & Iorigin[[i]] != "?") # This will be used to validate haploblocks based on the number of inf. alleles
        gaps = split(wildcardIndex, cumsum(c(1, diff(wildcardIndex) != 1)))
        if (length(gaps[[1]]) != 0 & length(wildcardIndex) != length(Iorigin[[i]])) {
          for (gap in gaps) {
            Iorigin[[i]][seq(gap[1], gap[length(gap)], 1)] = ifelse(gap[1] == 1,
                                                                    Iorigin[[i]][gap[length(gap)] + 1],
                                                                    ifelse(gap[length(gap)] == nSnp,
                                                                           Iorigin[[i]][gap[1] - 1],
                                                                           ifelse(Iorigin[[i]][gap[1] - 1] == Iorigin[[i]][gap[length(gap)] + 1],
                                                                                  Iorigin[[i]][gap[length(gap)] + 1],
                                                                                  "*")))
          }
        } else {
          Iorigin[[i]] = Iorigin[[i]] # If there is no haplotype to impute
        }
        # 6. Add the id haplotype to the haplotype origin matrix
        originMatrixTMP = rbind(originMatrixTMP, Iorigin[[i]])
        # 7.1 PHASE IMPUTATION
        # Depending on the selected imputation method, TH, MD and ME will be imputed
        wildcardIndex = which(Iphase[[i]] == "0")
        if (length(wildcardIndex) != 0) {
          if (imp == "imputeNot") {
            Iphase[[i]][wildcardIndex] = "0"
          } else if (imp == "imputeTHonly") {
            Iphase[[i]][wildcardIndex] = ifelse(noTripleHet[wildcardIndex] == FALSE,
                                                ifelse(Iorigin[[i]][wildcardIndex] == "P", 
                                                       Fphase[[i]][wildcardIndex],
                                                       ifelse(Iorigin[[i]][wildcardIndex] == "M",
                                                              Mphase[[i]][wildcardIndex],
                                                              "0")),
                                                "0")
          } else if (imp == "imputeAll") { 
            Iphase[[i]][wildcardIndex] = ifelse(Iorigin[[i]][wildcardIndex] == "P", 
                                                Fphase[[i]][wildcardIndex],
                                                ifelse(Iorigin[[i]][wildcardIndex] == "M",
                                                       Mphase[[i]][wildcardIndex],
                                                       "0"))
          }
        } else {
          Iphase[[i]] = Iphase[[i]]        
        }
        phaseMatrixTMP = rbind(phaseMatrixTMP, Iphase[[i]])
      }
      rowNames = append(rowNames, paste0(id, c("_P", "_M")))
      if (imp != "imputeNot" || cor != "correctNot") {
        # 7.2. PHASE CORRECTION FOR TRIPLE HETEROZYGOUS
        # Triple Het phases must be always heterozygous, despite imputation returning homozygous phases. We can also use the information of assigned alleles to assign 0s
        if (cor == "reImpute" || cor == "correctAll") {
          wildcardIndexP = which(phaseMatrixTMP[nrow(phaseMatrixTMP)-1, ] == 0)
          phaseMatrixTMP[nrow(phaseMatrixTMP)-1, wildcardIndexP] = ifelse(!noTripleHet[wildcardIndexP],
                                                                          ifelse(phaseMatrixTMP[nrow(phaseMatrixTMP), wildcardIndexP] == "1",
                                                                                 "2",
                                                                                 ifelse(phaseMatrixTMP[nrow(phaseMatrixTMP), wildcardIndexP] == "2",
                                                                                        "1",
                                                                                        "0")),
                                                                          "0")
          wildcardIndexM = which(phaseMatrixTMP[nrow(phaseMatrixTMP), ] == "0")
          phaseMatrixTMP[nrow(phaseMatrixTMP), wildcardIndexM] = ifelse(!noTripleHet[wildcardIndexM],
                                                                        ifelse(phaseMatrixTMP[nrow(phaseMatrixTMP)-1, wildcardIndexM] == "1",
                                                                               "0",
                                                                               ifelse(phaseMatrixTMP[nrow(phaseMatrixTMP)-1, wildcardIndexM] == "2",
                                                                                      "1",
                                                                                      "0")),
                                                                        "0") # For ME and MD, we will keep 0 if they lie within RE gaps to avoid phasing errors
        }
        if (cor == "correctFalseHom" || cor == "correctAll") {
          tripleHet = which(noTripleHet == FALSE)
          phaseMatrixTMPcopy = phaseMatrixTMP
          phaseMatrixTMP[nrow(phaseMatrixTMP)-1, tripleHet] = ifelse(phaseMatrixTMPcopy[nrow(phaseMatrixTMPcopy)-1, tripleHet] == phaseMatrixTMPcopy[nrow(phaseMatrixTMPcopy), tripleHet], # If pat and mat allele are hom in TH, 0
                                                                     "0",
                                                                     phaseMatrixTMP[nrow(phaseMatrixTMP)-1, tripleHet])
          phaseMatrixTMP[nrow(phaseMatrixTMP), tripleHet] = ifelse(phaseMatrixTMPcopy[nrow(phaseMatrixTMPcopy)-1, tripleHet] == phaseMatrixTMPcopy[nrow(phaseMatrixTMPcopy), tripleHet], # If pat and mat allele are hom in TH, 0
                                                                   "0",
                                                                   phaseMatrixTMP[nrow(phaseMatrixTMP), tripleHet])
          
        }
      }
    }
    rownames(originMatrixTMP) = rowNames ; rownames(phaseMatrixTMP) = rowNames ; rownames(originMatrix.noFilt.noImpTMP) = rowNames ; rownames(originMatrix.noFiltTMP) = rowNames
    originMatrix = rbind(originMatrix, originMatrixTMP)
    phaseMatrix = rbind(phaseMatrix, phaseMatrixTMP)
    originMatrix.noFilt.noImp = rbind(originMatrix.noFilt.noImp, originMatrix.noFilt.noImpTMP)
    originMatrix.noFilt = rbind(originMatrix.noFilt, originMatrix.noFiltTMP)
  }
  write.table(phaseMatrix, paste0(pop, "_", chr, ".phase"), sep = " ", row.names = TRUE, quote = FALSE, col.names = FALSE)
  write.table(originMatrix, paste0(pop, "_", chr, ".f.i.origin"), sep = " ", row.names = TRUE, quote = FALSE, col.names = FALSE)
  write.table(originMatrix.noFilt.noImp, paste0(pop, "_", chr, ".nf.ni.origin"), sep = " ", row.names = TRUE, quote = FALSE, col.names = FALSE)
  write.table(originMatrix.noFilt, paste0(pop, "_", chr, ".nf.i.origin"), sep = " ", row.names = TRUE, quote = FALSE, col.names = FALSE)
}
