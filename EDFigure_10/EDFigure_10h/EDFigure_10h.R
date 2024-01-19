
library(DiffBind)
library(rtracklayer)
library(GenomicRanges)
library(genomation)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(AnnotationDbi)
library(AnnotationHub)
library(ensembldb)
library(circlize)
library(profileplyr)
library(EnrichedHeatmap)
library(tidyverse)

# Astrocytes <- readBed("/Users/magalhae/Desktop/Julian_ChIPseq_NGN2/astrocyte/Astrocytes_hg38.bed", remove.unusual = T)
# Fetal_brain <- readBed("/Users/magalhae/Desktop/Julian_ChIPseq_NGN2/bed/Fetal_brain_hg38.bed", remove.unusual = T)
# ESC_neuron <- readBed("/Users/magalhae/Desktop/Julian_ChIPseq_NGN2/bed/ESC_neuron_hg38.bed", remove.unusual = T)
# 
# TAD = read.csv("/Users/magalhae/Desktop/Julian_ChIPseq_NGN2/bed/TAD_hg38.bed", sep = "\t")
# TE_package = read.csv("/Users/magalhae/Desktop/Julian_ChIPseq_NGN2/bed/TAD_hg38.bed", sep = "\t")
# SE_ele_package = read.csv("/Users/magalhae/Desktop/Julian_ChIPseq_NGN2/bed/SE_ele_package_hg38.bed", sep = "\t")
# SE_package = read.csv("/Users/magalhae/Desktop/Julian_ChIPseq_NGN2/bed/SE_package_hg38.bed", sep = "\t")


samples = read.csv('samples.csv')

dbObj <- dba(sampleSheet=samples) 

dbObj <- dba.blacklist(dbObj)


dba.save(dbObj, file='Chip_dbObj_PostCountCountNGN2', dir='.', pre='dba_', ext='RData', bRemoveAnalysis=FALSE, bRemoveBackground=FALSE, bCompress=FALSE)



dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE, bParallel = T)


dbObj <- dba.normalize(dbObj, normalize =DBA_NORM_LIB)


dbObj$config$RunParallel <- T


dba.save(dbObj, file='Chip_dbObj_NGN2_Norm', dir='.', pre='dba_', ext='RData', bRemoveAnalysis=FALSE, bRemoveBackground=FALSE, bCompress=FALSE)

#dbObj <- dba.load(file='Chip_dbObj_NGN2', dir='.', pre='dba_', ext='RData')


dba.plotPCA(dbObj,  attributes=c(DBA_FACTOR, DBA_CONDITION), label=DBA_REPLICATE)

