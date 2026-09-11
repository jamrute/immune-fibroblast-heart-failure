################################################################################
# CITE-seq: cell QC, RNA (SCTransform) and protein normalization
#
# Paper : Amrute, Luo et al., Nature 635, 423-433 (2024) | doi:10.1038/s41586-024-08008-5
#
# Purpose
#   Filters cells on gene/UMI counts and mitochondrial fraction, normalizes RNA
#   with SCTransform and runs PCA on RNA, CLR-normalized ADT and dsb protein
#   (excluding isotype controls and failed antibodies).
#
# Inputs
#   preprocessed.rds  (from 01_dsb_normalization.Rmd)
#
# Outputs
#   normalized.rds
#
# Run order
#   Upstream  : 01_dsb_normalization.Rmd
#   Downstream: 03_harmony_integration.Rmd
################################################################################

library(Seurat)
library(dplyr)

## ---- Paths (EDIT) ----
in_rds  <- "path/to/preprocessed.rds"   # EDIT: from 01_dsb_normalization.Rmd
out_rds <- "path/to/normalized.rds"     # EDIT

sample <- readRDS(in_rds)
sample <- UpdateSeuratObject(sample)

# QC
sample <- subset(sample,subset=nFeature_RNA>500&nFeature_RNA<6000&nCount_RNA>1000&nCount_RNA<25000&propmt<0.15) 

# Normalize RNA: SCTransform 
DefaultAssay(sample) <- 'RNA'
sample <- SCTransform(sample, vars.to.regress = c("propmt", "nCount_RNA"))
sample <- RunPCA(sample, verbose=FALSE, reduction.name="sctpca", npcs=100)

# Normalize Protein ADT and run PCA
DefaultAssay(sample) <- 'ADT'
VariableFeatures(sample) <- rownames(sample[["ADT"]])
sample <- NormalizeData(sample, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'adtpca')

# PCA on Protein DSB
# Exclude isotype controls and antibodies without staining (by row index)
prots = rownames(sample@assays$CITE@data)[-c(150,181,220:222,242:245)]
DefaultAssay(sample) <- 'CITE'
VariableFeatures(sample) <- prots
sample = ScaleData(sample, assay = 'CITE', verbose = FALSE)
sample = RunPCA(sample, reduction.name = 'pdsb', features = VariableFeatures(sample), verbose = FALSE)

saveRDS(sample, file = out_rds)