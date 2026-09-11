library(Seurat)
library(dplyr)
library(harmony)
#############################################
#Load scRNA files
#############################################

matrix <- Read10X(data.dir = "/data/Junedh/Amgen_ICM/Multiome/counts/Results/MA5/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "HF")
cells <- readRDS("/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/samplePostQC/MA5_cells")
s1 <- subset(scRNA, cells = cells)
s1$sample <- "MA5"

matrix <- Read10X(data.dir = "/data/Junedh/Amgen_ICM/Multiome/counts/Results/MA6/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "HF")
cells <- readRDS("/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/samplePostQC/MA6_cells")
s2 <- subset(scRNA, cells = cells)
s2$sample <- "MA6"

matrix <- Read10X(data.dir = "/data/Junedh/Amgen_ICM/Multiome/counts/Results/MA7/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "HF")
cells <- readRDS("/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/samplePostQC/MA7_cells")
s3 <- subset(scRNA, cells = cells)
s3$sample <- "MA7"

matrix <- Read10X(data.dir = "/data/Junedh/Amgen_ICM/Multiome/counts/Results/MA8/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "HF")
cells <- readRDS("/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/samplePostQC/MA8_cells")
s4 <- subset(scRNA, cells = cells)
s4$sample <- "MA8"

matrix <- Read10X(data.dir = "/data/Junedh/Amgen_ICM/Multiome/counts/Results/MA9/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "HF")
cells <- readRDS("/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/samplePostQC/MA9_cells")
s5 <- subset(scRNA, cells = cells)
s5$sample <- "MA9"

matrix <- Read10X(data.dir = "/data/Junedh/Amgen_ICM/Multiome/counts/Results/MA10/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "HF")
cells <- readRDS("/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/samplePostQC/MA10_cells")
s6 <- subset(scRNA, cells = cells)
s6$sample <- "MA10"

matrix <- Read10X(data.dir = "/data/Junedh/Amgen_ICM/Multiome/counts/Results/MA11/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "HF")
cells <- readRDS("/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/samplePostQC/MA11_cells")
s7 <- subset(scRNA, cells = cells)
s7$sample <- "MA11"

matrix <- Read10X(data.dir = "/data/Junedh/Amgen_ICM/Multiome/counts/Results/MA13/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "HF")
cells <- readRDS("/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/samplePostQC/MA13_cells")
s8 <- subset(scRNA, cells = cells)
s8$sample <- "MA13"

matrix <- Read10X(data.dir = "/data/Junedh/Amgen_ICM/Multiome/counts/Results/MA14/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "HF")
cells <- readRDS("/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/samplePostQC/MA14_cells")
s9 <- subset(scRNA, cells = cells)
s9$sample <- "MA14"

matrix <- Read10X(data.dir = "/data/Junedh/Amgen_ICM/Multiome/counts/batch3/MA19/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "HF")
cells <- readRDS("/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/samplePostQC/MA19_cells")
s10 <- subset(scRNA, cells = cells)
s10$sample <- "MA19"

matrix <- Read10X(data.dir = "/data/Junedh/Amgen_ICM/Multiome/counts/batch3/MA20/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "HF")
cells <- readRDS("/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/samplePostQC/MA20_cells")
s11 <- subset(scRNA, cells = cells)
s11$sample <- "MA20"

matrix <- Read10X(data.dir = "/data/Junedh/Amgen_ICM/Multiome/counts/batch3/MA22/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "HF")
cells <- readRDS("/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/samplePostQC/MA22_cells")
s12 <- subset(scRNA, cells = cells)
s12$sample <- "MA22"

matrix <- Read10X(data.dir = "/data/Junedh/Amgen_ICM/Multiome/counts/batch3/MA23/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "HF")
cells <- readRDS("/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/samplePostQC/MA23_cells")
s13 <- subset(scRNA, cells = cells)
s13$sample <- "MA23"

matrix <- Read10X(data.dir = "/data/Junedh/Amgen_ICM/Multiome/counts/batch3/MA24/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "HF")
cells <- readRDS("/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/samplePostQC/MA24_cells")
s14 <- subset(scRNA, cells = cells)
s14$sample <- "MA24"

matrix <- Read10X(data.dir = "/data/Junedh/Amgen_ICM/Multiome/counts/batch3/MA25/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "HF")
cells <- readRDS("/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/samplePostQC/MA25_cells")
s15 <- subset(scRNA, cells = cells)
s15$sample <- "MA25"

matrix <- Read10X(data.dir = "/data/Junedh/Amgen_ICM/Multiome/counts/batch3/MA26/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "HF")
cells <- readRDS("/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/samplePostQC/MA26_cells")
s16 <- subset(scRNA, cells = cells)
s16$sample <- "MA26"

matrix <- Read10X(data.dir = "/data/Junedh/Amgen_ICM/Multiome/counts/batch3/MA27/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "HF")
cells <- readRDS("/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/samplePostQC/MA27_cells")
s17 <- subset(scRNA, cells = cells)
s17$sample <- "MA27"

matrix <- Read10X(data.dir = "/data/Junedh/Amgen_ICM/Multiome/counts/batch3_additional_sequencing/MA28/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "HF")
cells <- readRDS("/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/samplePostQC/MA28_cells")
s18 <- subset(scRNA, cells = cells)
s18$sample <- "MA28"

matrix <- Read10X(data.dir = "/data/Junedh/Amgen_ICM/Multiome/counts/batch3_additional_sequencing/MA29/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "HF")
cells <- readRDS("/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/samplePostQC/MA29_cells")
s19 <- subset(scRNA, cells = cells)
s19$sample <- "MA29"

matrix <- Read10X(data.dir = "/data/Junedh/Amgen_ICM/Multiome/counts/batch3/MA30/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "HF")
cells <- readRDS("/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/samplePostQC/MA30_cells")
s20 <- subset(scRNA, cells = cells)
s20$sample <- "MA30"

matrix <- Read10X(data.dir = "/data/Junedh/Amgen_ICM/Multiome/counts/batch3/MA31/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "HF")
cells <- readRDS("/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/samplePostQC/MA31_cells")
s21 <- subset(scRNA, cells = cells)
s21$sample <- "MA31"

matrix <- Read10X(data.dir = "/data/Junedh/Amgen_ICM/Multiome/counts/batch3/MA32/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "HF")
cells <- readRDS("/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/samplePostQC/MA32_cells")
s22 <- subset(scRNA, cells = cells)
s22$sample <- "MA32"

matrix <- Read10X(data.dir = "/data/Junedh/Amgen_ICM/Multiome/counts/batch3/MA33/outs/filtered_feature_bc_matrix")	#Directory, not .h5 file
scRNA <- CreateSeuratObject(counts = matrix$`Gene Expression`, min.cells = 0, min.features = 0, project = "HF")
cells <- readRDS("/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/samplePostQC/MA33_cells")
s23 <- subset(scRNA, cells = cells)
s23$sample <- "MA33"

rm(matrix)
rm(scRNA)
rm(cells)

#Create merged Seurat object
merged <- merge(s1, y = c(s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21,s22,s23), 
                add.cell.ids = c("s1","s2","s3","s4","s5","s6","s7","s8","s9","s10",
                	"s11","s12","s13","s14","s15","s16","s17","s18","s19","s20",
                	"s21","s22","s23"), project = "HF")

merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^MT-")
merged <- subset(merged, subset = percent.mt < 5 & nCount_RNA < 15000)

saveRDS(merged, "/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/RNA/globalObjectConstruction_batch1_2_3/RNA_merged_postQC.rds")

DefaultAssay(merged) <- 'RNA'
merged <- SCTransform(merged, vars.to.regress = c("percent.mt", "nCount_RNA"))
merged <- RunPCA(merged, npcs=100, verbose=TRUE)
merged <- RunHarmony(merged, c("sample"), reduction = "pca", reduction.save = "harmony", assay.use = "SCT")
merged <- RunUMAP(merged, reduction = "harmony", dims = 1:50)
merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:50)
merged <- FindClusters(merged, graph.name = "SCT_snn", algorithm = 3, resolution = c(0.5,0.6,0.7,0.8,0.9,1.0), verbose = TRUE)

saveRDS(merged, "/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/RNA/globalObjectConstruction_batch1_2_3/integrated_preClean.rds")

Idents(merged) <- "SCT_snn_res.0.5"
DefaultAssay(merged) <- 'SCT'
rna.markers <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(rna.markers, file ="/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/RNA/globalObjectConstruction_batch1_2_3/DE_RNA_SCT_snn_res.0.5.csv", quote = FALSE)

Idents(merged) <- "SCT_snn_res.0.7"
DefaultAssay(merged) <- 'SCT'
rna.markers <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(rna.markers, file ="/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/RNA/globalObjectConstruction_batch1_2_3/DE_RNA_SCT_snn_res.0.7.csv", quote = FALSE)

Idents(merged) <- "SCT_snn_res.1"
DefaultAssay(merged) <- 'SCT'
rna.markers <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(rna.markers, file ="/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/RNA/globalObjectConstruction_batch1_2_3/DE_RNA_SCT_snn_res.1.0.csv", quote = FALSE)







