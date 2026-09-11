library(ArchR)
addArchRGenome("hg38")
library(Seurat)
library(dplyr)

setwd(paste("/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/", sep=""))

scRNA <- readRDS("/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/RNA/globalObjectConstruction_batch1_2_3/integrated/myocardium_multiome.rds")
proj1 <- loadArchRProject(path = "/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/Save-proj1")

# Approach 1: Create peak2gene bed file
p2g <- getPeak2GeneLinks(
    ArchRProj = proj1,
    corCutOff = 0.3,
    resolution = 1,
    returnLoops = FALSE
)

p2gChr <- (metadata(p2g)$peakSet %>% seqnames(.))[p2g$idxATAC]
p2gStart <- (metadata(p2g)$peakSet %>% start(.))[p2g$idxATAC]
p2gEnd <- (metadata(p2g)$peakSet %>% end(.))[p2g$idxATAC]
p2gGene <- (metadata(p2g)$geneSet)$name[p2g$idxRNA]

p2gBed <- data.frame(p2gChr, p2gStart, p2gEnd, p2gGene)
write.table(p2gBed, "./Save-proj1/ArchR_output/p2g.bed", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)



# ##########################################################  GPC Analysis  ############################################################
# #GSM/GEX correlation
# corr_vals <- correlateMatrices(proj1, useMatrix1 = "GeneScoreMatrix", useMatrix2 = "GeneExpressionMatrix",reducedDims = "HAR_Combined") #replace with GeneExpressionMatrix for multiome
# corr_vals <- as.data.frame(corr_vals)

# #remove insignificant correlations
# corr_vals <- corr_vals[!is.na(corr_vals$cor),]
# corr_vals <- corr_vals[corr_vals$padj < 0.05,]

# #get peak 2 gene links
# p2g <- getPeak2GeneLinks(proj1, returnLoops = FALSE)
# x <- readRDS("/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/Save-proj1/Peak2GeneLinks/seRNA-Group-KNN.rds")
# p2g$idxRNA <- x@rowRanges$name[p2g$idxRNA]

# gpc <- as.data.frame(table(p2g$idxRNA))
# rownames(gpc) <- gpc$Var1

# gpc <- merge(gpc,corr_vals[c('GeneScoreMatrix_name','cor')], by.x='Var1', by.y='GeneScoreMatrix_name')

# #ggscatter(x='Freq',y='cor',data = gpc)

# write.csv(gpc,"./Save-proj1/ArchR_output/gpc_pilot_table.csv", quote = FALSE)
# write.csv(corr_vals,"./Save-proj1/ArchR_output/gex_gsm_corr_pilot_table.csv", quote = FALSE)