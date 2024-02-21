#!/usr/bin/env Rscript

#==================================================================================
#title: 0.haplomagic.R
#description: controls the execution of the pipeline
#author: jmontero
#email: jose.a.montero-tena@ab.uni-giessen.de
#date: 2023-02-21
#version: 1.0.0
#notes:
#==================================================================================

suppressMessages(suppressWarnings(if(!require("fread")){ install.packages("fread") }))

args = commandArgs(trailingOnly=TRUE)

# Define global variables
pops = unlist(strsplit(args[1], split = " "))
chrs = unlist(strsplit(args[2], split = " "))
min = as.numeric(args[3])
imp = as.character(args[4]) # 'imputeNot', 'imputeTHonly' or 'imputeAll'
cor = as.character(args[5]) # if 'imp' != 'correctNot', 'reImpute', 'correctFalseHom', 'correctAll'
thr = as.numeric(args[6])

# Source functions
source("1.phase.founders_function.R")
source("2.phase.rest_function.R")
source("3.haploblock_function.R")
source("4.detect.recombination_function.R")
source("5.statistics_function.R")

# Loop pipeline
for (pop in pops) {
  for (chr in chrs) {
    PhaseFounders(pop, chr)
    PhaseRest(pop, chr, min, imp, cor)
    AssignFounderHaploblocks(pop, chr)
    DetectRecombination(pop, chr, thr)
    GetHaplomagicStatistics(pop, chr)
  }
}
