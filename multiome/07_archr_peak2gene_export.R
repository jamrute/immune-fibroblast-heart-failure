################################################################################
# Multiome: export peak-to-gene links (BED) and gene-peak connectivity (GPC)
#
# Paper : Amrute, Luo et al., Nature 635, 423-433 (2024) | doi:10.1038/s41586-024-08008-5
# Part  : Single-nucleus multiome (paired RNA + ATAC)
#
# Purpose
#   Writes significant peak-to-gene links (r > 0.3) as a BED file and,
#   optionally, computes the number of linked peaks per gene together with the
#   correlation between gene activity scores and gene expression.
#
# Inputs  : ArchR project (Save-proj1) with peak-to-gene links
# Outputs : Save-proj1/ArchR_output/p2g.bed
#           (optional) gpc_pilot_table.csv, gex_gsm_corr_pilot_table.csv
#
# Run order
#   Upstream  : 06_archr_peak2gene.R
#   Downstream: none
################################################################################

library(ArchR)
addArchRGenome("hg38")
library(Seurat)
library(dplyr)

## ---- Paths and options (EDIT) ----
analysis_dir <- "path/to/analysis_batch_1_2_3"   # EDIT: folder containing Save-proj1/
run_gpc      <- FALSE  # TRUE to also compute the gene-peak connectivity tables
setwd(analysis_dir)

proj1 <- loadArchRProject(path = "Save-proj1")

## ---- Peak-to-gene links as BED (chr, start, end, gene) ----
p2g <- getPeak2GeneLinks(ArchRProj = proj1, corCutOff = 0.3, resolution = 1, returnLoops = FALSE)

p2gChr   <- (metadata(p2g)$peakSet %>% seqnames(.))[p2g$idxATAC]
p2gStart <- (metadata(p2g)$peakSet %>% start(.))[p2g$idxATAC]
p2gEnd   <- (metadata(p2g)$peakSet %>% end(.))[p2g$idxATAC]
p2gGene  <- (metadata(p2g)$geneSet)$name[p2g$idxRNA]

p2gBed <- data.frame(p2gChr, p2gStart, p2gEnd, p2gGene)
write.table(p2gBed, "./Save-proj1/ArchR_output/p2g.bed", sep = "\t", col.names = FALSE,
            row.names = FALSE, quote = FALSE)

## ---- Gene-peak connectivity (GPC) ----
if (run_gpc) {
  # Correlation between gene activity (GeneScoreMatrix) and expression, significant only
  corr_vals <- correlateMatrices(proj1, useMatrix1 = "GeneScoreMatrix", useMatrix2 = "GeneExpressionMatrix",
                                 reducedDims = "HAR_Combined")
  corr_vals <- as.data.frame(corr_vals)
  corr_vals <- corr_vals[!is.na(corr_vals$cor), ]
  corr_vals <- corr_vals[corr_vals$padj < 0.05, ]

  # Number of linked peaks per gene
  p2g <- getPeak2GeneLinks(proj1, returnLoops = FALSE)
  x <- readRDS("./Save-proj1/Peak2GeneLinks/seRNA-Group-KNN.rds")
  p2g$idxRNA <- x@rowRanges$name[p2g$idxRNA]
  gpc <- as.data.frame(table(p2g$idxRNA))
  rownames(gpc) <- gpc$Var1
  gpc <- merge(gpc, corr_vals[c('GeneScoreMatrix_name', 'cor')], by.x = 'Var1', by.y = 'GeneScoreMatrix_name')

  write.csv(gpc, "./Save-proj1/ArchR_output/gpc_pilot_table.csv", quote = FALSE)
  write.csv(corr_vals, "./Save-proj1/ArchR_output/gex_gsm_corr_pilot_table.csv", quote = FALSE)
}
