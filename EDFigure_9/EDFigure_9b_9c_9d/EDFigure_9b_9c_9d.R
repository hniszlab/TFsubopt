---
  title: "R Notebook"
output: html_notebook
---
  
  This is an [R Markdown](http://rmarkdown.rstudio.com) Notebook. When you execute code within the notebook, the results appear beneath the code. 

Try executing this chunk by clicking the *Run* button within the chunk or by placing your cursor inside it and pressing *Cmd+Shift+Enter*. 


library(DiffBind)
library(rtracklayer)
library(GenomicRanges)
library(genomation)
library(profileplyr)
library(EnrichedHeatmap)
library(circlize)

library(tidyverse)
library(ChIPQC)
#library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(ChIPseeker)
library(clusterProfiler)
library(AnnotationDbi)
library(AnnotationHub)
library(ensembldb)
library(org.Mm.eg.db)
library(Rsamtools)
library(Signac)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

setwd("/Users/magalhae/Desktop/Julian_ChIP_CEBPA/")


### This files were supplyed by  https://doi.org/10.7554/eLife.65381

my.file=system.file("BLaER_SE.bed",package="genomation")
refseq = readBed("BLaER_SE.bed",track.line=FALSE,remove.unusual=FALSE)


### This files were supplyed by  https://doi.org/10.7554/eLife.65381

BLaER_S_gr =  readBed("BLaER_SE_cut.bed",track.line=FALSE,remove.unusual=FALSE)
enhancers_gr =  readBed("enhancers_sorted_cut.bed",track.line=FALSE,remove.unusual=FALSE)
iMac_SE_gr =  readBed("iMac_SE_cut.bed",track.line=FALSE,remove.unusual=FALSE)
promoters_gr =  readBed("promoters_sorted_cut.bed",track.line=FALSE,remove.unusual=FALSE)
common_gr = readBed("Common_Sites_cut.bed",track.line=FALSE,remove.unusual=FALSE)

### This files were supplyed by  https://doi.org/10.7554/eLife.65381

BLaER_S_obj =  read.table("BLaER_SE.bed", header = FALSE)
enhancers_obj =  read.table("enhancers_sorted.bed", header = FALSE)
iMac_SE_obj =  read.table("iMac_SE.bed", header = FALSE)
promoters_obj =  read.table("promoters_sorted.bed", header = FALSE)
common_obj = read.table("Common_Sites.bed", header = FALSE)

# split and convert per region
BLaER_S <- 
  lapply(split(BLaER_S_obj, BLaER_S_obj$V4), function(i){
    GRanges(seqnames = i$V1,
            ranges = IRanges(start = i$V2,
                             end = i$V3,
                             names = i$V4))
  })

# split and convert per region
enhancers_l <- 
  lapply(split(enhancers_obj, enhancers_obj$V4), function(i){
    GRanges(seqnames = i$V1,
            ranges = IRanges(start = i$V2,
                             end = i$V3,
                             names = i$V4))
  })
# split and convert per region
iMac_SE <- 
  lapply(split(iMac_SE_obj, iMac_SE_obj$V4), function(i){
    GRanges(seqnames = i$V1,
            ranges = IRanges(start = i$V2,
                             end = i$V3,
                             names = i$V4))
  })
# split and convert per region
promoters_l <- 
  lapply(split(promoters_obj, promoters_obj$V4), function(i){
    GRanges(seqnames = i$V1,
            ranges = IRanges(start = i$V2,
                             end = i$V3,
                             names = i$V4))
  })

# split and convert per region
common_l <- 
  lapply(split(common_obj, common_obj$V4), function(i){
    GRanges(seqnames = i$V1,
            ranges = IRanges(start = i$V2,
                             end = i$V3,
                             names = i$V4))
  })









dbObj <- dba(sampleSheet=samples)



dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE)


dba.plotPCA(dbObj,  attributes=c(DBA_FACTOR, DBA_CONDITION), label=DBA_ID)


p <- plot(dbObj)


dbObj <- dba.contrast(dbObj, categories=DBA_FACTOR, minMembers = 2)


dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)


t <- dba.show(dbObj, bContrasts=T)
write_csv(t, "Chipseq_CEBPa_contrasts.csv")



dba.show(dbObj, bContrasts=T)




dba.plotPCA(dbObj, contrast=1, method=DBA_DESEQ2, attributes=DBA_FACTOR, label=DBA_ID)



dba.plotVenn(dbObj,contrast=2,method=DBA_ALL_METHODS)
dba.plotVenn(dbObj,contrast=5,method=DBA_ALL_METHODS)
dba.plotVenn(dbObj,contrast=1,method=DBA_ALL_METHODS)
dba.plotVenn(dbObj,contrast=6,method=DBA_ALL_METHODS)


dba.plotBox(dbObj, contrast=1, attributes=c(DBA_FACTOR, DBA_CONDITION), notch=FALSE)
dba.plotBox(dbObj, contrast=2, attributes=c(DBA_FACTOR, DBA_CONDITION), notch=FALSE)
dba.plotBox(dbObj, contrast=5, attributes=c(DBA_FACTOR, DBA_CONDITION), notch=FALSE)
dba.plotBox(dbObj, contrast=6, attributes=c(DBA_FACTOR, DBA_CONDITION), notch=FALSE)


res_deseq1 <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)
out_t1 <- annotatePeak(res_deseq1,  TxDb=txdb, tssRegion=c(-2000, 2000))
out_df1 <- data.frame(out_t1@anno)

res_deseq2 <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 2, th=1)
out_t2 <- annotatePeak(res_deseq2,  TxDb=txdb, tssRegion=c(-2000, 2000))
out_df2 <- data.frame(out_t2@anno)

res_deseq5 <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 5, th=1)
out_t5 <- annotatePeak(res_deseq5,  TxDb=txdb, tssRegion=c(-2000, 2000))
out_df5 <- data.frame(out_t5@anno)

res_deseq6 <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 6, th=1)
out_t6 <- annotatePeak(res_deseq6,  TxDb=txdb, tssRegion=c(-2000, 2000))
out_df6 <- data.frame(out_t6@anno)

entrez = NULL
# Get the entrez IDs

entrez <- out_df1$geneId
entrez <- out_df2$geneId
entrez <- out_df5$geneId
entrez <- out_df6$geneId



entrez <- unique(entrez)

# Return the gene symbol for the set of Entrez IDs
annotations_edb <- AnnotationDbi::select(ahEdb,
                                         keys = entrez,
                                         columns = c("GENENAME"),
                                         keytype = "ENTREZID")

# Change IDs to character type to merge
annotations_edb$ENTREZID <- as.character(annotations_edb$ENTREZID)

out1 <- out_df1 %>% left_join(annotations_edb, by=c("geneId"="ENTREZID"))
out2 <- out_df2 %>% left_join(annotations_edb, by=c("geneId"="ENTREZID"))
out5 <- out_df5 %>% left_join(annotations_edb, by=c("geneId"="ENTREZID"))
out6 <- out_df6 %>% left_join(annotations_edb, by=c("geneId"="ENTREZID"))

write_csv(out1, "Diffbind_results_ChIP_CEBPa_CEBPa_IS15_24vsCEBPa_IS15_48.csv")
write_csv(out2, "Diffbind_results_ChIP_CEBPa_CEBPa_IS15_24vsCEBPa_WT_24.csv")
write_csv(out5, "Diffbind_results_ChIP_CEBPa_CEBPa_IS15_48vsCEBPa_WT_48.csv")
write_csv(out6, "Diffbind_results_ChIP_CEBPa_CEBPa_WT_24vsCEBPa_WT_48.csv")





#profiles <- dba.plotProfile(dbObj,samples= dbObj[["samples"]][["SampleID"]], merge=DBA_REPLICATE)
#profiles_p <- dba.plotProfile(dbObj,samples= dbObj[["samples"]][["SampleID"]], merge=DBA_REPLICATE, sites=promoters_gr)

profiles_p2 <- dba.plotProfile(dbObj,samples= dbObj[["samples"]][["SampleID"]], merge=DBA_REPLICATE, sites=promoters_gr, nOfWindows=27, distanceAround=100 , style="percentOfRegion")
profiles_e <- dba.plotProfile(dbObj,samples= dbObj[["samples"]][["SampleID"]], merge=DBA_REPLICATE, sites=enhancers_gr)

profiles_e2 <- dba.plotProfile(dbObj,samples= dbObj[["samples"]][["SampleID"]], merge=DBA_REPLICATE, sites=enhancers_gr, nOfWindows=27, distanceAround=100 , style="percentOfRegion")


profiles_bSE2 <- dba.plotProfile(dbObj,samples= dbObj[["samples"]][["SampleID"]], merge=DBA_REPLICATE, sites=BLaER_S_gr, nOfWindows=27, distanceAround=100 , style="percentOfRegion")
profiles_mSE2 <- dba.plotProfile(dbObj,samples= dbObj[["samples"]][["SampleID"]], merge=DBA_REPLICATE, sites=iMac_SE_gr, nOfWindows=27, distanceAround=100 , style="percentOfRegion")
#profiles_comon <- dba.plotProfile(dbObj,samples= dbObj[["samples"]][["SampleID"]], merge=DBA_REPLICATE, sites=common_gr)


pdf("profiles_CEBPa.pdf")
dba.plotProfile(profiles, use_raster = T, raster_device = "CairoTIFF", raster_quality = 10, matrices_color = col_list)
dev.off()

pdf("profiles_promoter_CEBPa.pdf")
dba.plotProfile(profiles_p, use_raster = T, raster_device = "CairoTIFF", raster_quality = 10, matrices_color = col_list)
dev.off()

pdf("profiles_promoter_CEBPa_bined.pdf")
dba.plotProfile(profiles_p2, use_raster = T, raster_device = "CairoTIFF", raster_quality = 10, matrices_color = col_list)
dev.off()

pdf("profiles_enhancer_CEBPa.pdf")
dba.plotProfile(profiles_e, use_raster = T, raster_device = "CairoTIFF", raster_quality = 10, matrices_color = col_list)
dev.off()

pdf("profiles_bcell_SE_CEBPa_binded.pdf")
dba.plotProfile(profiles_bSE2, use_raster = T, raster_device = "CairoTIFF", raster_quality = 10, matrices_color = col_list)
dev.off()

pdf("profiles_macro_SE_CEBPa_binded.pdf")
dba.plotProfile(profiles_mSE2, use_raster = T, raster_device = "CairoTIFF", raster_quality = 10, matrices_color = col_list)
dev.off()

pdf("profiles_enhancer_CEBPa_binded.pdf")
dba.plotProfile(profiles_e2, use_raster = T, raster_device = "CairoTIFF", raster_quality = 10, matrices_color = col_list)
dev.off()

pdf("profiles_ComWTIS15_CEBPa.pdf")
dba.plotProfile(profiles_comon, use_raster = T, raster_device = "CairoTIFF", raster_quality = 10, matrices_color = col_list)
dev.off()


profilestest  <- dba.plotProfile(dbObj, sites=2)
profilestest5  <- dba.plotProfile(dbObj, sites=5)

repObj <- dba.report(dbObj, contrast=c(2,5), bDB=TRUE, bGain=T, bLoss=T, bAll=T)
repObj2 <- dba.report(dbObj, contrast=c(1,2,5,6), bDB=TRUE, bGain=T, bLoss=T)
repObj3 <- dba.report(dbObj, contrast=c(2,5), bDB=TRUE, bGain=T, bLoss=T, bAll=F)

profiles_dif <- dba.plotProfile(dbObj,samples= dbObj[["samples"]][["SampleID"]], merge=DBA_REPLICATE, sites=repObj)
profiles_dif2 <- dba.plotProfile(dbObj,samples= dbObj[["samples"]][["SampleID"]], merge=DBA_REPLICATE, sites=repObj2)
profiles_dif3 <- dba.plotProfile(dbObj,samples= dbObj[["samples"]][["SampleID"]], merge=DBA_REPLICATE, sites=repObj3)




dba.show(dbObj, bContrast=T)


#repList <- GRangesList(Gain=repObj[repObj$peaks[2]],Loss=repObj[repObj$peaks[5]])

repList <- GRangesList(Gain24=repObj[["peaks"]][[5]],Gain48=repObj[["peaks"]][[2]],Loss24=repObj[["peaks"]][[3]],Loss48=repObj[["peaks"]][[6]])


profiles_v2 <- dba.plotProfile(dbObj, sites=repList)



dba.plotProfile(profiles_v2, use_raster = T, raster_device = "CairoTIFF", raster_quality = 2, matrices_color = col_list)


dba.plotProfile(profiles, use_raster = T, raster_device = "CairoTIFF", raster_quality = 2)





overlaps <- dba.plotVenn(dbObj,contrast=c(2,5))



profiles_v <- dba.plotProfile(dbObj, sites=overlaps)


generateEnrichedHeatmap(profiles_bSE2)


pdf("profiles_dif.pdf")
dba.plotProfile(profiles_dif, use_raster = T, raster_device = "CairoTIFF", raster_quality = 10, matrices_color = col_list)
dev.off()


pdf("profiles_Comon.pdf")
dba.plotProfile(profiles_comon, use_raster = T, raster_device = "CairoTIFF", raster_quality = 10, matrices_color = col_list)
dev.off()


pdf("profiles_dif2.pdf")
dba.plotProfile(profiles_v, use_raster = T, raster_device = "CairoTIFF", raster_quality = 10, matrices_color = col_list)
dev.off()




col_fun_24IS15 <- colorRamp2(c(0, 30), c("white", "#652d90"))
col_fun_48IS15 <- colorRamp2(c(0, 30), c("white", "#652d90"))
col_fun_24WT <- colorRamp2(c(0, 30), c("white", "#f05756"))
col_fun_48WT <- colorRamp2(c(0, 30), c("white", "#f05756"))

col_list = c(col_fun_24IS15, col_fun_48IS15, col_fun_24WT, col_fun_48WT)
names(col_list) <- c("CEBPa_IS15_24_24", "CEBPa_IS15_48_48", "CEBPa_WT_24_24", "CEBPa_WT_48_48" )


generateEnrichedHeatmap(profiles_v, color_by_sample_group = "sample_labels", all_color_scales_equal = TRUE, matrices_color = col_list)


overlapReg <- dba.overlap(dbObj,dbObj$masks$All,mode=DBA_OLAP_ALL,bCorOnly=F)



profiles_dif3 <- dba.plotProfile(dbObj,samples= list(CEBPa_IS15_24=dbObj$mask$CEBPa_IS15_24, CEBPa_IS15_48=dbObj$mask$CEBPa_IS15_48,
                                                     CEBPa_WT_24=dbObj$mask$CEBPa_WT_24, CEBPa_WT_48=dbObj$mask$CEBPa_WT_48), merge=DBA_REPLICATE, sites= t)



profiles_dif2 <- dba.plotProfile(dbObj,samples= list(CEBPa_IS15_24=dbObj$mask$CEBPa_IS15_24, CEBPa_IS15_48=dbObj$mask$CEBPa_IS15_48, 
                                                     CEBPa_WT_24=dbObj$mask$CEBPa_WT_24, CEBPa_WT_48=dbObj$mask$CEBPa_WT_48), merge=DBA_REPLICATE, sites=repObj2)


dba.plotProfile(profiles_dif, use_raster = T, raster_device = "CairoTIFF", raster_quality = 2)




output_path <- file.path(tempdir(),"profiles_dif.MAT.gz")
output_path2 <- file.path(tempdir(),"profiles_dif2.MAT.gz")
output_pathc <- file.path(tempdir(),"profiles_comon.MAT.gz")

#export_deepToolsMat(profiles_dif, con = output_path)
export_deepToolsMat(profiles_dif2, con = output_path2)
export_deepToolsMat(profiles_comon, con = output_pathc)

EH_matdf2_5_3 <- convertToEnrichedHeatmapMat(profiles_dif3)
EH_matdf2_5 <- convertToEnrichedHeatmapMat(profiles_dif)
EH_matdf1_2_5_6 <- convertToEnrichedHeatmapMat(profiles_dif2)
EH_matdfcom <- convertToEnrichedHeatmapMat(profiles_comon)


heatmap <- generateEnrichedHeatmap(profiles_dif3)


heatmap_Test <- generateEnrichedHeatmap(profiles_diftest)



class <- repObj[["class"]]


col_fun_24IS15 <- colorRamp2(c(0, 30), c("white", "#652d90"))
col_fun_48IS15 <- colorRamp2(c(0, 30), c("white", "#652d90"))
col_fun_24WT <- colorRamp2(c(0, 30), c("white", "#f05756"))
col_fun_48WT <- colorRamp2(c(0, 30), c("white", "#f05756"))

col_list = c(col_fun_24IS15, col_fun_48IS15, col_fun_24WT, col_fun_48WT)
names(col_list) <- c("CEBPa_IS15_24_24", "CEBPa_IS15_48_48", "CEBPa_WT_24_24", "CEBPa_WT_48_48" )


heatmap2 <- generateEnrichedHeatmap(profiles_dif, color_by_sample_group = "sample_labels", all_color_scales_equal = TRUE, matrices_color = col_list)



pdf("plotProfile_ISvsWT_gainloss.pdf")
dba.plotProfile(profiles_dif, use_raster = T, raster_device = "CairoTIFF", raster_quality = 10)
dev.off()
pdf("plotProfile_ISvsWT_48vs24_gainloss.pdf")
dba.plotProfile(profiles_dif2, use_raster = T, raster_device = "CairoTIFF", raster_quality = 10)
dev.off()




dba.plotProfile(profilestest, use_raster = T, raster_device = "CairoTIFF", raster_quality = 8)
dba.plotProfile(profiles_dif, use_raster = T, raster_device = "CairoTIFF", raster_quality = 8)


dba.plotProfile(profiles_dif, use_raster = T, raster_device = "CairoTIFF", raster_quality = 8)
dba.plotProfile(profiles_dif2, use_raster = T, raster_device = "CairoTIFF", raster_quality = 8)



dba.plotProfile(profiles_dif3, use_raster = T, raster_device = "CairoTIFF", raster_quality = 2)



