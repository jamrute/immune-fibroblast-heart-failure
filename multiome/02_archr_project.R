library(ArchR)
addArchRGenome("hg38")
library(Seurat)
library(dplyr)

################ QC and Generate the ArchR project after filtering and add gene expression matrix from paired scRNAseq data

#Get input fragment files for each sample
arrow_loc <- "/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/samplePostQC/"

samples <- c("MA5","MA6","MA7","MA8","MA9","MA10","MA11","MA13","MA14","MA19","MA20","MA22","MA23","MA24","MA25",
             "MA26","MA27","MA28","MA29","MA30","MA31","MA32","MA33")

ArrowFiles <- c()

for (s in samples) {
  ArrowFiles <- c(ArrowFiles, paste(arrow_loc, s, "/", s, ".arrow", sep=""))
}

proj1 <- ArchRProject(ArrowFiles, copyArrows = TRUE)

scRNA <- readRDS("/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/RNA/globalObjectConstruction_batch1_2_3/integrated/myocardium_multiome.rds")
col <- readRDS("/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/RNA/globalObjectConstruction_batch1_2_3/integrated/myocardium_multiome_cellList")

#Isolate cells used in Seurat
proj1 <- subsetCells(ArchRProj = proj1, cellNames = col)

# GEX Matrix
gene_matrix_files <- c()

b1 <- c("MA5","MA6","MA7","MA8","MA9","MA10","MA11","MA13","MA14")
b2 <- c("MA19","MA20","MA22","MA23","MA24","MA25","MA26","MA27","MA30","MA31","MA32","MA33")
b3 <- c("MA28","MA29")


# RNA input matrix
for (s in samples) {
  if (s %in% b1) {
  rna_file <- c(paste("/data/Junedh/Amgen_ICM/Multiome/counts/Results/", s, "/outs/filtered_feature_bc_matrix.h5", sep=""))
  } else if (s %in% b2) {
  rna_file <- c(paste("/data/Junedh/Amgen_ICM/Multiome/counts/batch3/", s, "/outs/filtered_feature_bc_matrix.h5", sep=""))
  } else if (s %in% b3) {
  rna_file <- c(paste("/data/Junedh/Amgen_ICM/Multiome/counts/batch3_additional_sequencing/", s, "/outs/filtered_feature_bc_matrix.h5", sep=""))
  }

  gene_matrix_files <- c(gene_matrix_files, rna_file)
}

seRNA <- import10xFeatureMatrix(input = gene_matrix_files, names = samples)

seRNAcombined<-cbind(assay(seRNA[[1]]), assay(seRNA[[2]]), assay(seRNA[[3]]), assay(seRNA[[4]]), assay(seRNA[[5]]), assay(seRNA[[6]]), assay(seRNA[[7]]),
                     assay(seRNA[[8]]), assay(seRNA[[9]]), assay(seRNA[[10]]), assay(seRNA[[11]]), assay(seRNA[[12]]), assay(seRNA[[13]]), assay(seRNA[[14]]),
                     assay(seRNA[[14]]), assay(seRNA[[15]]), assay(seRNA[[16]]), assay(seRNA[[17]]), assay(seRNA[[18]]), assay(seRNA[[19]]), assay(seRNA[[20]]),
                     assay(seRNA[[21]]), assay(seRNA[[22]]), assay(seRNA[[23]]))

seRNA2<-SummarizedExperiment(assays=list(counts=seRNAcombined), rowRanges= rowRanges(seRNA[[1]]))

proj1 <- addGeneExpressionMatrix(
  input = proj1,
  seRNA = seRNA2,
  chromSizes = getChromSizes(proj1),
  excludeChr = c("chrM", "chrY"),
  scaleTo = 10000,
  verbose = TRUE,
  threads = getArchRThreads(),
  parallelParam = NULL,
  force = TRUE,
  logFile = createLogFile("addGeneExpressionMatrix")
)

# Transfer cluster and identities from scRNA
mappedCellType <- as.character(scRNA$predicted.celltype)
HFetiology <- as.character(scRNA$HF.etiology)
condition <- as.character(scRNA$condition)

proj1$mappedCellType <- mappedCellType
proj1$HFetiology <- HFetiology
proj1$condition <- condition

proj1 <- saveArchRProject(ArchRProj = proj1, outputDirectory = "/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/Save-proj1", load = TRUE)




