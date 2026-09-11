################################################################################
# Shared helper functions used across the analysis scripts
#
# Paper : Amrute, Luo et al., Nature 635, 423-433 (2024) | doi:10.1038/s41586-024-08008-5
#
# Usage : source(file.path(repo_dir, "R", "utils.R"))
#         where `repo_dir` points to the root of this repository.
#
# These functions replace blocks that were previously copy-pasted many times.
# gene_set_zscore(), annotate_clusters() and plot_composition() were tested to
# reproduce the original code exactly.
################################################################################

suppressPackageStartupMessages(library(ggplot2))

# ---- Scoring, annotation and plotting ----------------------------------------

#' Gene-set z-score
#'
#' Scales each gene across cells (z-score), sets genes with zero variance to 0,
#' and averages the z-scores of the matched genes for every cell. Gene names are
#' matched case-insensitively, so one list works for human and mouse symbols.
#'
#' @param expdata Gene x cell matrix, e.g. GetAssayData(obj) for the default assay.
#' @param genes Character vector of gene symbols.
#' @return Named numeric vector (one score per cell).
gene_set_zscore <- function(expdata, genes) {
  zz <- which(tolower(rownames(expdata)) %in% tolower(genes))
  if (length(zz) == 0) stop("None of the genes were found in the expression matrix.")
  geneExp <- as.matrix(expdata[zz, ])
  geneExp <- t(scale(t(geneExp)))
  geneExp[is.nan(geneExp)] <- 0
  colSums(geneExp) / length(zz)
}

#' Map cluster IDs to annotation labels
#'
#' @param clusters Vector of cluster IDs (e.g. obj$SCT_snn_res.0.5); names are kept.
#' @param labels Named character vector, e.g. c("0" = "Fibroblast", "1" = "Myeloid").
#' @return Character vector of labels. Stops if a cluster has no label.
annotate_clusters <- function(clusters, labels) {
  ids <- as.character(clusters)
  missing <- setdiff(unique(ids), names(labels))
  if (length(missing) > 0) stop("No label for cluster(s): ", paste(missing, collapse = ", "))
  stats::setNames(unname(labels[ids]), names(clusters))
}

#' Stacked proportion bar plot (cell composition per group)
#'
#' @param meta Cell metadata (e.g. obj@meta.data).
#' @param x,fill Column names for the x axis and the fill.
#' @param colors Fill colours passed to scale_fill_manual().
plot_composition <- function(meta, x, fill, colors) {
  ggplot(meta, aes(x = !!rlang::sym(x), fill = !!rlang::sym(fill))) +
    geom_bar(position = "fill") + theme_linedraw() +
    theme(axis.text.x = element_text(angle = 90)) +
    scale_fill_manual(values = colors) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
}

# ---- Loading data ---------------------------------------------------------------

#' Read a sample sheet (CSV) with every column as character
read_sample_sheet <- function(path) {
  utils::read.csv(path, colClasses = "character", check.names = FALSE)
}

#' Load and SCTransform-normalize Visium sections listed in a sample sheet
#'
#' @param sheet Data frame with a `file` column (h5seurat file name) plus metadata
#'   columns (e.g. sample, sample_id, region, days_postMI, disease) that are added
#'   to every spot of that section.
#' @param data_dir Folder containing the h5seurat files.
#' @return Named list of Seurat objects (names = `sample` column).
load_visium_samples <- function(sheet, data_dir) {
  objs <- lapply(seq_len(nrow(sheet)), function(i) {
    obj <- SeuratDisk::LoadH5Seurat(file.path(data_dir, sheet$file[i]))
    for (col in setdiff(names(sheet), "file")) obj[[col]] <- sheet[[col]][i]
    Seurat::SCTransform(obj, assay = "RNA", verbose = FALSE)
  })
  stats::setNames(objs, sheet$sample)
}

#' Convert an h5ad file to h5seurat (next to it) and load it
load_h5ad_as_seurat <- function(h5ad) {
  SeuratDisk::Convert(h5ad, dest = "h5seurat", overwrite = TRUE)
  SeuratDisk::LoadH5Seurat(sub("\\.h5ad$", ".h5seurat", h5ad))
}

#' Paths to the Cell Ranger ARC filtered matrices for each multiome sample
#'
#' @param sheet Data frame with columns `sample` and `batch_dir`.
#' @param counts_dir Root folder holding the Cell Ranger ARC runs.
#' @param h5 TRUE for the .h5 file, FALSE for the matrix directory.
multiome_matrix_paths <- function(sheet, counts_dir, h5 = FALSE) {
  leaf <- if (h5) "filtered_feature_bc_matrix.h5" else "filtered_feature_bc_matrix"
  file.path(counts_dir, sheet$batch_dir, sheet$sample, "outs", leaf)
}

#' Load multiome gene-expression counts, keeping only cells that passed ATAC QC
#'
#' @param sheet Data frame with columns `sample` and `batch_dir`.
#' @param counts_dir Root folder holding the Cell Ranger ARC runs.
#' @param cells_dir Folder with one `<sample>_cells` RDS file (barcodes) per sample.
#' @return Named list of Seurat objects with a `sample` metadata column.
load_multiome_rna <- function(sheet, counts_dir, cells_dir, project = "HF") {
  paths <- multiome_matrix_paths(sheet, counts_dir)
  objs <- lapply(seq_len(nrow(sheet)), function(i) {
    counts <- Seurat::Read10X(data.dir = paths[i])
    obj <- Seurat::CreateSeuratObject(counts = counts$`Gene Expression`, min.cells = 0,
                                      min.features = 0, project = project)
    obj <- subset(obj, cells = readRDS(file.path(cells_dir, paste0(sheet$sample[i], "_cells"))))
    obj$sample <- sheet$sample[i]
    obj
  })
  stats::setNames(objs, sheet$sample)
}

# ---- Differential testing ---------------------------------------------------------

#' Positive cluster markers at several clustering resolutions, one CSV each
#'
#' @param obj Seurat object with clustering columns named <prefix><resolution>.
#' @param resolutions Named character vector: names are the resolution labels used
#'   in the metadata column, values are the labels used in the output file name.
write_cluster_markers <- function(obj, resolutions, out_dir, prefix = "SCT_snn_res.",
                                  assay = "SCT") {
  for (res in names(resolutions)) {
    Seurat::Idents(obj) <- paste0(prefix, res)
    Seurat::DefaultAssay(obj) <- assay
    markers <- Seurat::FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
    utils::write.csv(markers, file = file.path(out_dir, paste0("DE_RNA_SCT_snn_res.", resolutions[[res]], ".csv")),
                     quote = FALSE)
  }
}

#' Differential peak accessibility between two groups (ArchR)
#'
#' Runs a Wilcoxon test on the PeakMatrix (TSS enrichment and fragment-count
#' bias matched), saves MA and volcano plots and the result object.
#'
#' @param proj ArchRProject.
#' @param group_by cellColData column holding the groups (e.g. "cellType_condition").
#' @param use_group,bgd_group Foreground and background group names.
#' @param out_dir Folder for the RDS result.
#' @return The getMarkerFeatures() result (invisibly).
archr_differential_peaks <- function(proj, group_by, use_group, bgd_group,
                                     out_dir = "./Save-proj1/ArchR_output",
                                     cutoff = "FDR <= 0.1 & abs(Log2FC) >= 0.5") {
  marker_test <- ArchR::getMarkerFeatures(
    ArchRProj = proj, useMatrix = "PeakMatrix", groupBy = group_by,
    testMethod = "wilcoxon", bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = use_group, bgdGroups = bgd_group, maxCells = 5000
  )
  pma <- ArchR::markerPlot(seMarker = marker_test, name = use_group, cutOff = cutoff, plotAs = "MA")
  pv  <- ArchR::markerPlot(seMarker = marker_test, name = use_group, cutOff = cutoff, plotAs = "Volcano")
  label <- paste0(use_group, "-vs-", bgd_group)
  ArchR::plotPDF(pma, pv, name = paste0(label, "-Markers-MA-Volcano"), width = 5, height = 5,
                 ArchRProj = proj, addDOC = FALSE)
  saveRDS(marker_test, file.path(out_dir, paste0(label, ".rds")))
  invisible(marker_test)
}

#' Row-normalized confusion-matrix heatmap between two cell labelings (ArchR)
plot_confusion_heatmap <- function(labels_row, labels_col) {
  cM <- ArchR::confusionMatrix(paste0(labels_row), paste0(labels_col))
  cM <- cM / Matrix::rowSums(cM)
  pheatmap::pheatmap(mat = as.matrix(cM), color = ArchR::paletteContinuous("whiteBlue"),
                     border_color = "black")
}
