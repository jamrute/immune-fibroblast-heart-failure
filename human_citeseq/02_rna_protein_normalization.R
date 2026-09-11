library(Seurat)
library(dplyr)

sample <- readRDS("/data/Junedh/Amgen_ICM/Human/preprocessed.rds")
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
prots = rownames(sample@assays$CITE@data)[-c(150,181,220:222,242:245)]
DefaultAssay(sample) <- 'CITE'
VariableFeatures(sample) <- prots
sample = ScaleData(sample, assay = 'CITE', verbose = FALSE)
sample = RunPCA(sample, reduction.name = 'pdsb', features = VariableFeatures(sample), verbose = FALSE)

saveRDS(sample, file = "/data/Junedh/Amgen_ICM/Human/normalized.rds")