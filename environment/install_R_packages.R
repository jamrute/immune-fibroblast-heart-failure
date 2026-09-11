################################################################################
# Install the R packages used in this repository.
# Paper : Amrute, Luo et al., Nature 635, 423-433 (2024) | doi:10.1038/s41586-024-08008-5
#
# Usage : Rscript environment/install_R_packages.R
# Notes : Scripts use Seurat v4-style assay access (obj@assays$RNA@counts,
#         slot = "data"); with Seurat v5 set
#         options(Seurat.object.assay.version = "v3") before creating objects.
################################################################################

options(repos = c(CRAN = "https://cloud.r-project.org"))
for (p in c("BiocManager", "remotes")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

cran <- c("Seurat", "sctransform", "harmony", "dsb", "dplyr", "tidyr", "tibble", "tidyverse", "purrr",
          "magrittr", "data.table", "Matrix", "reshape2", "stringr", "ggplot2", "ggpubr", "ggsci",
          "patchwork", "cowplot", "pheatmap", "viridis", "RColorBrewer", "MetBrewer", "scales",
          "colorspace", "png")
bioc <- c("SingleCellExperiment", "S4Vectors", "scater", "DESeq2", "edgeR", "apeglm", "clusterProfiler",
          "DOSE", "enrichplot", "ReactomePA", "biomaRt", "progeny", "dorothea", "viper", "Nebulosa",
          "BiocGenerics")
github <- c("GreenleafLab/ArchR",          # multiome ATAC analysis; paletteDiscrete()/paletteContinuous()
            "mojaveazure/seurat-disk",     # SeuratDisk: h5Seurat <-> h5ad
            "cvarrichio/Matrix.utils")     # aggregate.Matrix() for pseudobulk (archived on CRAN)

install.packages(setdiff(cran, rownames(installed.packages())))
BiocManager::install(setdiff(bioc, rownames(installed.packages())), update = FALSE)
for (repo in github) remotes::install_github(repo, upgrade = "never")
# ArchR also needs: ArchR::installExtraPackages(); MACS2 on the PATH for peak calling

message("Done. Record versions with sessionInfo() for reproducibility.")
