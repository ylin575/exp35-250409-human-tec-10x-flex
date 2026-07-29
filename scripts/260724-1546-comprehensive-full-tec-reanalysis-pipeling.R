# =============================================================================
# HUMAN TEC 10x FLEX -- COMPREHENSIVE RE-ANALYSIS PIPELINE
#
# PURPOSE: rebuild the analysis from raw count matrices with every QC and
# artifact check we've established, then test under 4 conditions whether the
# "cluster 8" phenomenon (cTEC/mTEC-I marker co-expression) persists:
#
#   Condition 1: all samples,  no Harmony   <- reproduces original analysis
#   Condition 2: all samples,  + Harmony    <- does integration remove it?
#   Condition 3: drop HT70,    no Harmony   <- is it HT70-dependent?
#   Condition 4: drop HT70,    + Harmony    <- both controls together
#
# WHY HT70 SPECIFICALLY:
#   - only male donor (drives Y-gene signal in cluster 8 sub-regions)
#   - 70.7% of near_mTEC_I, 71.4% of near_Nurse (vs 29.1% dataset baseline)
#   - SHORTER FLEX FIXATION: 10 h vs 16 h for HT67 and HT71
#     (10x recommends overnight fixation; 10 h is a real protocol deviation)
#   Three independent reasons HT70 could be an outlier -- leave-one-out tests
#   all three at once without needing to disentangle them.
#
# ORIGINAL PARAMETERS PRESERVED (from 260506-no-integration-code-clean-up.R):
#   min.features = 200 at object creation
#   percent.mt <= mean + 3*SD  (per sample)
#   nFeature_RNA >= 500 & <= mean + 3*SD  (per sample)
#   FindNeighbors(dims = 1:18); FindClusters(resolution = 0.5); RunUMAP(dims = 1:18)
#
# DELIBERATE CHANGES FROM ORIGINAL (each justified inline):
#   min.cells = 0 per sample, gene filtering applied ONCE after merge
#   add.cell.ids = sample names (interpretable barcodes)
#   scDblFinder doublet detection added
#   cell cycle scoring added
# =============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(scCustomize)
library(patchwork)
library(harmony)
library(scDblFinder)
library(SingleCellExperiment)

DATE_PREFIX <- format(Sys.time(), "%y%m%d-%H%M")
set.seed(42)   # scDblFinder + UMAP involve randomness

OUTDIR <- "results"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
out <- function(x) file.path(OUTDIR, paste0(DATE_PREFIX, "-", x))

# =============================================================================
# GLOBAL CONFIG
# =============================================================================
H5_DIR         <- "rawdata/h5-files"
MIN_FEATURES   <- 200     # at object creation (original)
QC_MIN_FEATURE <- 500     # QC floor (original)
SD_THRESHOLD   <- 3       # 3*SD upper bounds (original)
POOLED_MIN_CELLS <- 3     # gene filter, applied AFTER merge (see rationale)
PCA_DIMS       <- 1:18    # original
CLUSTER_RES    <- 0.5     # original
REMOVE_DOUBLETS <- TRUE   # set FALSE to keep doublets flagged but retained

SAMPLE_IDS <- c("ht67-cd205neg", "ht67-cd205pos", "ht70", "ht71")

# Donor metadata (from sample-info.xlsx)
SAMPLE_META <- data.frame(
  sample.id      = SAMPLE_IDS,
  donor          = c("HT67", "HT67", "HT70", "HT71"),
  age_days       = c(90, 90, 365, 1460),
  sex            = c("Female", "Female", "Male", "Female"),
  institution    = c("Mount Sinai", "Mount Sinai", "Mount Sinai", "Northwell"),
  fixation_hours = c(16, 16, 10, 16),
  sort_strategy  = c("CD205neg", "CD205pos", "AllTEC", "AllTEC"),
  stringsAsFactors = FALSE
)

# Expected cell counts from 10x-flex-cell-counts.xlsx (for sanity checking)
EXPECTED_CELLS <- c("ht67-cd205neg" = 8280, "ht67-cd205pos" = 12953,
                    "ht70" = 13246, "ht71" = 11073)

# -----------------------------------------------------------------------------
# Marker panels (from the manuscript + prior analyses in this project)
# -----------------------------------------------------------------------------
cTEC_markers    <- c("PSMB11", "PRSS16", "LY75", "TBATA")
mTEC_I_markers  <- c("KRT14", "CCL19", "KRT15", "KRT19", "FN1", "CLU")
mTEC_II_markers <- c("AIRE", "FEZF2")
BMC_markers     <- c("COL4A5", "COL4A6", "COL8A1", "ITGB4")
Nurse_markers   <- c("TBATA","PRSS16","PTPRC", "CD3E", "RAG1", "CD4")   # thymocyte-engulfment signature

# Cluster-8 signatures for tracking across runs. Cluster NUMBERS are NOT stable
# between runs -- we must identify the population by its transcriptional
# signature instead. Two signatures used:
#   (a) pooled cluster-8 markers (from FindMarkers vs all)
#   (b) near_cTEC-specific markers -- the sub-region the reviewer asked about
C8_POOLED_SIG <- c("LRRC55", "POU2AF1", "CHST1", "FOXQ1", "RELT", "FATE1",
                   "FOXS1", "PROM2", "TLCD1", "FXYD2", "TRAF1", "STAR",
                   "DAPK2", "SORCS1", "CALML3")
C8_NEARCTEC_SIG <- c("CHST1", "FOXS1", "KCNE5", "FOXQ1", "CCL5", "SEZ6L2",
                     "NPHS1", "PROM2", "LTBP2", "TLCD1", "CRYAB", "CALML3",
                     "RELT", "ENPP6", "DAPK2")

# Sex-linked genes for donor-confound checking
SEX_GENES <- c("XIST", "RPS4Y1", "DDX3Y", "UTY", "KDM5D", "EIF1AY", "USP9Y")


# =============================================================================
# SANITY CHECK HELPER
# =============================================================================
sanity <- function(label, ..., fatal = TRUE) {
  cat("\n--- SANITY CHECK:", label, "---\n")
  checks <- list(...)
  all_ok <- TRUE
  for (nm in names(checks)) {
    ok <- isTRUE(checks[[nm]])
    cat(sprintf("  [%s] %s\n", ifelse(ok, "PASS", "FAIL"), nm))
    if (!ok) all_ok <- FALSE
  }
  if (!all_ok && fatal) stop("Sanity check failed at: ", label)
  invisible(all_ok)
}


# =============================================================================
# PART A -- LOAD, QC, AND BUILD THE BASE OBJECT (done once, shared by all runs)
# =============================================================================

# -----------------------------------------------------------------------------
# A1. Load count matrices
# -----------------------------------------------------------------------------
cat("\n", strrep("=", 70), "\nPART A1: Loading count matrices\n", strrep("=", 70), "\n")

raw_mats <- lapply(SAMPLE_IDS, function(s) {
  path <- file.path(H5_DIR, s, "sample_filtered_feature_bc_matrix.h5")
  cat("Reading:", path, "\n")
  Read10X_h5(path, use.names = TRUE, unique.features = TRUE)
})
names(raw_mats) <- SAMPLE_IDS

# SANITY: same gene panel across all samples? (10x Flex uses one probe set;
# differing gene counts here would indicate a probe-panel version mismatch)
gene_counts <- sapply(raw_mats, nrow)
cell_counts_raw <- sapply(raw_mats, ncol)
cat("\nGenes per sample:\n"); print(gene_counts)
cat("\nCells per sample (raw h5):\n"); print(cell_counts_raw)
cat("\nExpected cells (from cell-count spreadsheet):\n"); print(EXPECTED_CELLS)

sanity("A1 -- raw matrices loaded",
       "all samples loaded"            = length(raw_mats) == 4,
       "identical gene panel across samples" = length(unique(gene_counts)) == 1,
       "cell counts within 5% of expected"   = all(
         abs(cell_counts_raw[names(EXPECTED_CELLS)] - EXPECTED_CELLS) /
           EXPECTED_CELLS < 0.05)
)

# -----------------------------------------------------------------------------
# A2. Create Seurat objects
# -----------------------------------------------------------------------------
# CHANGE FROM ORIGINAL: min.cells = 0 instead of 3.
# RATIONALE: per-sample min.cells filtering removes a DIFFERENT set of genes
# from each sample (a gene detected in 2 cells in one sample but 50 in another
# survives in one and not the other). On merge, the union is taken and the
# missing sample's cells are zero-filled -- creating artificial sample-
# correlated variance that can leak into HVG selection and PCA. Filtering once
# on the pooled object (step A6) applies the same bar to every sample.
cat("\n", strrep("=", 70), "\nPART A2: Creating Seurat objects\n", strrep("=", 70), "\n")

seu_list <- lapply(SAMPLE_IDS, function(s) {
  obj <- CreateSeuratObject(
    counts       = raw_mats[[s]],
    project      = "htec-10xflex",
    min.cells    = 0,               # <- pooled filtering instead (A6)
    min.features = MIN_FEATURES
  )
  obj$sample.id  <- s
  obj$percent.mt <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj
})
names(seu_list) <- SAMPLE_IDS

cells_after_create <- sapply(seu_list, ncol)
cat("\nCells after min.features =", MIN_FEATURES, ":\n"); print(cells_after_create)
cat("Cells dropped:\n"); print(cell_counts_raw - cells_after_create)

sanity("A2 -- Seurat objects created",
       "no sample lost all cells"     = all(cells_after_create > 0),
       "no sample lost >50% of cells" = all(cells_after_create / cell_counts_raw > 0.5),
       "percent.mt computed"          = all(sapply(seu_list, function(o) "percent.mt" %in% colnames(o@meta.data))),
       "percent.mt in valid range"    = all(sapply(seu_list, function(o) all(o$percent.mt >= 0 & o$percent.mt <= 100)))
)

# -----------------------------------------------------------------------------
# A3. Per-sample QC filtering (thresholds identical to original script)
# -----------------------------------------------------------------------------
cat("\n", strrep("=", 70), "\nPART A3: Per-sample QC filtering\n", strrep("=", 70), "\n")

qc_thresholds <- list()
seu_list_qc <- lapply(SAMPLE_IDS, function(s) {
  obj <- seu_list[[s]]
  upper_mt      <- mean(obj$percent.mt)   + SD_THRESHOLD * sd(obj$percent.mt)
  upper_feature <- mean(obj$nFeature_RNA) + SD_THRESHOLD * sd(obj$nFeature_RNA)
  qc_thresholds[[s]] <<- c(upper_mt = upper_mt, upper_feature = upper_feature)
  
  obj <- subset(obj, subset = percent.mt <= upper_mt)
  obj <- subset(obj, subset = nFeature_RNA >= QC_MIN_FEATURE &
                  nFeature_RNA <= upper_feature)
  obj
})
names(seu_list_qc) <- SAMPLE_IDS

cat("\nQC thresholds applied (mean + 3*SD, per sample):\n")
print(do.call(rbind, qc_thresholds))

cells_after_qc <- sapply(seu_list_qc, ncol)
cat("\nCells after QC:\n"); print(cells_after_qc)
cat("Percent retained:\n")
print(round(100 * cells_after_qc / cells_after_create, 1))

sanity("A3 -- QC filtering",
       "all samples retain cells"           = all(cells_after_qc > 0),
       "all samples retain >60% of cells"   = all(cells_after_qc / cells_after_create > 0.60),
       "no sample disproportionately lost"  = (max(cells_after_qc / cells_after_create) -
                                                 min(cells_after_qc / cells_after_create)) < 0.25
)

# -----------------------------------------------------------------------------
# A4. Doublet detection -- PER SAMPLE (doublets only form within a GEM well)
# -----------------------------------------------------------------------------
cat("\n", strrep("=", 70), "\nPART A4: Doublet detection (scDblFinder)\n", strrep("=", 70), "\n")

seu_list_dbl <- lapply(SAMPLE_IDS, function(s) {
  cat("\nRunning scDblFinder on:", s, "\n")
  obj <- seu_list_qc[[s]]
  sce <- as.SingleCellExperiment(obj, assay = "RNA")
  sce <- scDblFinder(sce, verbose = TRUE)
  # coerce explicitly -- SCE returns a factor, which breaks dplyr downstream
  obj$scDblFinder_class <- as.character(sce$scDblFinder.class)
  obj$scDblFinder_score <- as.numeric(sce$scDblFinder.score)
  obj
})
names(seu_list_dbl) <- SAMPLE_IDS

doublet_rates <- sapply(seu_list_dbl, function(o)
  round(100 * mean(o$scDblFinder_class == "doublet"), 2))
cat("\nDoublet rate (%) per sample:\n"); print(doublet_rates)

sanity("A4 -- doublet detection",
       "all samples classified"        = all(sapply(seu_list_dbl, function(o)
         !any(is.na(o$scDblFinder_class)))),
       "doublet rates plausible (<25%)" = all(doublet_rates < 25),
       "doublet rates non-zero"         = all(doublet_rates > 0)
)

if (REMOVE_DOUBLETS) {
  seu_list_dbl <- lapply(seu_list_dbl, function(o)
    subset(o, subset = scDblFinder_class == "singlet"))
  cat("\nDoublets REMOVED. Cells remaining:\n")
  print(sapply(seu_list_dbl, ncol))
} else {
  cat("\nDoublets FLAGGED but retained (REMOVE_DOUBLETS = FALSE).\n")
}

# -----------------------------------------------------------------------------
# A5. Merge
# -----------------------------------------------------------------------------
# CHANGE FROM ORIGINAL: add.cell.ids supplied. The original merge() had none,
# so Seurat appended positional suffixes (_1.._4) -- which made mapping cells
# back to samples error-prone (this bit us during the SoupX barcode matching).
cat("\n", strrep("=", 70), "\nPART A5: Merging samples\n", strrep("=", 70), "\n")

samples <- merge(
  x = seu_list_dbl[[1]],
  y = seu_list_dbl[-1],
  add.cell.ids = SAMPLE_IDS
)
samples <- JoinLayers(samples)

expected_total <- sum(sapply(seu_list_dbl, ncol))
cat("\nMerged object:", ncol(samples), "cells x", nrow(samples), "genes\n")
cat("Cells per sample after merge:\n"); print(table(samples$sample.id))

sanity("A5 -- merge",
       "total cells = sum of parts"   = ncol(samples) == expected_total,
       "all 4 samples present"        = length(unique(samples$sample.id)) == 4,
       "barcodes unique"              = !any(duplicated(colnames(samples))),
       "layers joined (single counts layer)" = "counts" %in% Layers(samples, assay = "RNA")
)

# -----------------------------------------------------------------------------
# A6. Pooled gene filtering (the min.cells step, applied symmetrically)
# -----------------------------------------------------------------------------
cat("\n", strrep("=", 70), "\nPART A6: Pooled gene filtering\n", strrep("=", 70), "\n")

counts_mat <- GetAssayData(samples, assay = "RNA", layer = "counts")
genes_keep <- rownames(counts_mat)[Matrix::rowSums(counts_mat > 0) >= POOLED_MIN_CELLS]
cat("Genes before:", nrow(samples), " | after (>=", POOLED_MIN_CELLS,
    "cells):", length(genes_keep), "\n")

samples <- subset(samples, features = genes_keep)

sanity("A6 -- pooled gene filter",
       "genes retained"                = nrow(samples) == length(genes_keep),
       "cell count unchanged"          = ncol(samples) == expected_total,
       "key marker genes survived"     = all(c(cTEC_markers, mTEC_I_markers,
                                               mTEC_II_markers) %in% rownames(samples))
)

# -----------------------------------------------------------------------------
# A7. Attach donor metadata + cell cycle scoring
# -----------------------------------------------------------------------------
cat("\n", strrep("=", 70), "\nPART A7: Metadata + cell cycle\n", strrep("=", 70), "\n")

idx <- match(samples$sample.id, SAMPLE_META$sample.id)
samples$donor          <- SAMPLE_META$donor[idx]
samples$age_days       <- SAMPLE_META$age_days[idx]
samples$sex            <- SAMPLE_META$sex[idx]
samples$institution    <- SAMPLE_META$institution[idx]
samples$fixation_hours <- SAMPLE_META$fixation_hours[idx]
samples$sort_strategy  <- SAMPLE_META$sort_strategy[idx]

samples <- NormalizeData(samples, verbose = TRUE)   # needed for CellCycleScoring
samples <- CellCycleScoring(
  samples,
  s.features   = cc.genes.updated.2019$s.genes,
  g2m.features = cc.genes.updated.2019$g2m.genes,
  set.ident    = FALSE, verbose = TRUE
)

cat("\nMetadata cross-check:\n")
print(unique(samples@meta.data[, c("sample.id", "donor", "sex",
                                   "age_days", "fixation_hours")]))
cat("\nCell cycle phase distribution:\n")
print(table(samples$Phase, samples$sample.id))

sanity("A7 -- metadata + cell cycle",
       "no NA in donor metadata"  = !any(is.na(samples$donor)) && !any(is.na(samples$sex)),
       "Phase assigned to all"    = !any(is.na(samples$Phase)),
       "exactly one male donor"   = length(unique(samples$donor[samples$sex == "Male"])) == 1
)

# Sex-gene sanity: confirm the male donor is the one showing Y-gene expression
sex_genes_present <- intersect(SEX_GENES, rownames(samples))
if (length(sex_genes_present) > 0) {
  p_sex <- DotPlot_scCustom(samples, features = sex_genes_present,
                            group.by = "sample.id") + RotatedAxis() +
    labs(title = "Sex-linked genes by sample (ground-truth check)")
  ggsave(out("A7-sex-gene-check.pdf"), p_sex, width = 7, height = 4)
}

saveRDS(samples, out("base-object-QCd-3SD.rds"))
cat("\nBase QC'd object saved. Cells:", ncol(samples), " Genes:", nrow(samples), "\n")


# =============================================================================
# PART B -- PIPELINE FUNCTION (runs one condition end-to-end)
# =============================================================================

run_pipeline <- function(base_obj, include_ht70, use_harmony, label) {
  
  cat("\n\n", strrep("#", 70), "\n")
  cat("# CONDITION:", label, "\n")
  cat("#   include HT70:", include_ht70, " | Harmony:", use_harmony, "\n")
  cat(strrep("#", 70), "\n")
  
  obj <- base_obj
  
  # --- B1. Optional HT70 removal ---
  if (!include_ht70) {
    n_before <- ncol(obj)
    obj <- subset(obj, subset = sample.id != "ht70")
    cat("\nHT70 removed:", n_before, "->", ncol(obj), "cells\n")
    sanity(paste(label, "-- HT70 removal"),
           "ht70 absent"         = !("ht70" %in% unique(obj$sample.id)),
           "3 samples remain"    = length(unique(obj$sample.id)) == 3,
           "cells actually dropped" = ncol(obj) < n_before
    )
  }
  
  # --- B2. Standard preprocessing ---
  obj <- NormalizeData(obj, verbose = TRUE)
  obj <- FindVariableFeatures(obj, verbose = TRUE)
  obj <- ScaleData(obj, verbose = TRUE)
  obj <- RunPCA(obj, verbose = TRUE)
  
  sanity(paste(label, "-- preprocessing"),
         "HVGs selected"      = length(VariableFeatures(obj)) > 0,
         "PCA computed"       = "pca" %in% Reductions(obj),
         "enough PCs for dims" = ncol(Embeddings(obj, "pca")) >= max(PCA_DIMS)
  )
  
  # --- B3. Optional Harmony integration ---
  reduction_use <- "pca"
  if (use_harmony) {
    obj <- RunHarmony(obj, group.by.vars = "sample.id", reduction = "pca",
                      dims.use = PCA_DIMS, reduction.save = "harmony",
                      verbose = TRUE)
    reduction_use <- "harmony"
    sanity(paste(label, "-- Harmony"),
           "harmony reduction created" = "harmony" %in% Reductions(obj),
           "dims preserved"            = ncol(Embeddings(obj, "harmony")) >= max(PCA_DIMS)
    )
  }
  
  # --- B4. Clustering + UMAP (original parameters) ---
  obj <- FindNeighbors(obj, reduction = reduction_use, dims = PCA_DIMS, verbose = TRUE)
  obj <- FindClusters(obj, resolution = CLUSTER_RES, verbose = TRUE)
  obj <- RunUMAP(obj, reduction = reduction_use, dims = PCA_DIMS, verbose = TRUE)
  
  n_clusters <- length(unique(obj$seurat_clusters))
  cat("\nClusters found:", n_clusters, "\n")
  print(table(obj$seurat_clusters))
  
  sanity(paste(label, "-- clustering"),
         "clusters formed"           = n_clusters > 1,
         "cluster count plausible"   = n_clusters >= 10 && n_clusters <= 50,
         "UMAP computed"             = "umap" %in% Reductions(obj),
         "every cell assigned"       = !any(is.na(obj$seurat_clusters))
  )
  
  # --- B5. Module scores (identity panels + cluster-8 signatures) ---
  score_module <- function(o, genes, name) {
    g <- intersect(genes, rownames(o))
    if (length(g) == 0) { o[[name]] <- 0; return(o) }
    o <- AddModuleScore(o, features = list(g), name = name, verbose = TRUE)
    o[[name]] <- o[[paste0(name, "1")]]
    o[[paste0(name, "1")]] <- NULL
    o
  }
  obj <- score_module(obj, cTEC_markers,     "cTEC_score")
  obj <- score_module(obj, mTEC_I_markers,   "mTEC_I_score")
  obj <- score_module(obj, mTEC_II_markers,  "mTEC_II_score")
  obj <- score_module(obj, BMC_markers,      "BMC_score")
  obj <- score_module(obj, Nurse_markers,    "Nurse_score")
  obj <- score_module(obj, C8_POOLED_SIG,    "C8_pooled_score")
  obj <- score_module(obj, C8_NEARCTEC_SIG,  "C8_nearcTEC_score")
  
  # --- B6. THE KEY TEST: does a cluster-8-like population exist? ---
  # Cluster NUMBERS are not comparable across runs, so identify by signature.
  #   (a) which cluster scores highest on the C8 signatures?
  #   (b) which cluster has the highest cTEC/mTEC-I double-positive rate?
  counts_now <- GetAssayData(obj, assay = "RNA", layer = "counts")
  cTEC_avail   <- intersect(cTEC_markers,   rownames(counts_now))
  mTEC_I_avail <- intersect(mTEC_I_markers, rownames(counts_now))
  
  cTEC_detected_genes   <- Matrix::colSums(counts_now[cTEC_avail, , drop = FALSE] > 0)
  mTEC_I_detected_genes <- Matrix::colSums(counts_now[mTEC_I_avail, , drop = FALSE] > 0)
  
  obj$double_pos_strict <- (cTEC_detected_genes >= ceiling(length(cTEC_avail) / 2)) &
    (mTEC_I_detected_genes >= ceiling(length(mTEC_I_avail) / 2))
  
  cluster_summary <- obj@meta.data %>%
    group_by(seurat_clusters) %>%
    summarise(
      n                  = n(),
      pct_double_pos     = round(100 * mean(double_pos_strict), 1),
      mean_C8_pooled     = round(mean(C8_pooled_score), 3),
      mean_C8_nearcTEC   = round(mean(C8_nearcTEC_score), 3),
      mean_cTEC          = round(mean(cTEC_score), 3),
      mean_mTEC_I        = round(mean(mTEC_I_score), 3),
      pct_ht70           = round(100 * mean(sample.id == "ht70"), 1),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_C8_nearcTEC))
  
  cat("\n=== Cluster summary, ranked by near_cTEC C8 signature ===\n")
  print(as.data.frame(cluster_summary))
  write.csv(cluster_summary, out(paste0(label, "-cluster-summary.csv")),
            row.names = FALSE)
  
  top_c8 <- cluster_summary$seurat_clusters[1]
  cat("\n>>> Top C8-signature cluster:", as.character(top_c8),
      "| n =", cluster_summary$n[1],
      "| double-positive =", cluster_summary$pct_double_pos[1], "%\n")
  
  # --- B7. Plots ---
  p1 <- DimPlot_scCustom(obj, label = TRUE) + labs(title = paste0(label, ": clusters"))
  p2 <- DimPlot_scCustom(obj, group.by = "sample.id") + labs(title = "by sample")
  p3 <- FeaturePlot_scCustom(obj, features = "C8_nearcTEC_score",
                             colors_use = viridis_light_high, order = TRUE,
                             na_cutoff = NA) + labs(title = "C8 near_cTEC signature")
  p4 <- FeaturePlot_scCustom(obj, features = "double_pos_strict",
                             colors_use = viridis_light_high, order = TRUE,
                             na_cutoff = NA) + labs(title = "cTEC+/mTEC-I+ double-positive")
  ggsave(out(paste0(label, "-overview.pdf")), (p1 | p2) / (p3 | p4),
         width = 14, height = 11)
  
  saveRDS(obj, out(paste0(label, "-object.rds")))
  
  list(
    label            = label,
    object           = obj,
    n_cells          = ncol(obj),
    n_clusters       = n_clusters,
    cluster_summary  = cluster_summary,
    top_c8_cluster   = as.character(top_c8),
    top_c8_n         = cluster_summary$n[1],
    top_c8_dblpos    = cluster_summary$pct_double_pos[1],
    top_c8_sig       = cluster_summary$mean_C8_nearcTEC[1],
    overall_dblpos   = round(100 * mean(obj$double_pos_strict), 1)
  )
}


DimPlot_scCustom(res4$object, cells.highlight = WhichCells(res4$object, idents = "0"))
DimPlot_scCustom(res4$object, cells.highlight = WhichCells(res4$object, idents = "1"))
DimPlot_scCustom(res4$object, cells.highlight = WhichCells(res4$object, idents = "2"))
DimPlot_scCustom(res4$object, cells.highlight = WhichCells(res4$object, idents = "3"))
DimPlot_scCustom(res4$object, cells.highlight = WhichCells(res4$object, idents = "4"))
DimPlot_scCustom(res4$object, cells.highlight = WhichCells(res4$object, idents = "5"))
DimPlot_scCustom(res4$object, cells.highlight = WhichCells(res4$object, idents = "6"))
DimPlot_scCustom(res4$object, cells.highlight = WhichCells(res4$object, idents = "7"))
DimPlot_scCustom(res4$object, cells.highlight = WhichCells(res4$object, idents = "8"))
DimPlot_scCustom(res4$object, cells.highlight = WhichCells(res4$object, idents = "9"))
DimPlot_scCustom(res4$object, cells.highlight = WhichCells(res4$object, idents = "10"))
DimPlot_scCustom(res4$object, cells.highlight = WhichCells(res4$object, idents = "11"))
DimPlot_scCustom(res4$object, cells.highlight = WhichCells(res4$object, idents = "12"))
DimPlot_scCustom(res4$object, cells.highlight = WhichCells(res4$object, idents = "13"))
DimPlot_scCustom(res4$object, cells.highlight = WhichCells(res4$object, idents = "14"))
DimPlot_scCustom(res4$object, cells.highlight = WhichCells(res4$object, idents = "15"))
DimPlot_scCustom(res4$object, cells.highlight = WhichCells(res4$object, idents = "16"))
DimPlot_scCustom(res4$object, cells.highlight = WhichCells(res4$object, idents = "17"))
DimPlot_scCustom(res4$object, cells.highlight = WhichCells(res4$object, idents = "18"))
DimPlot_scCustom(res4$object, cells.highlight = WhichCells(res4$object, idents = "19"))
DimPlot_scCustom(res4$object, cells.highlight = WhichCells(res4$object, idents = "20"))
DimPlot_scCustom(res4$object, cells.highlight = WhichCells(res4$object, idents = "21"))

# =============================================================================
# PART C -- RUN ALL FOUR CONDITIONS
# =============================================================================
res1 <- run_pipeline(samples, include_ht70 = TRUE,  use_harmony = FALSE, label = "C1-all-noHarmony")
res2 <- run_pipeline(samples, include_ht70 = TRUE,  use_harmony = TRUE,  label = "C2-all-Harmony")
res3 <- run_pipeline(samples, include_ht70 = FALSE, use_harmony = FALSE, label = "C3-noHT70-noHarmony")
res4 <- run_pipeline(samples, include_ht70 = FALSE, use_harmony = TRUE,  label = "C4-noHT70-Harmony")

# =============================================================================
# PART C (extension) -- CONDITIONS 5 & 6: isolate the CD205-sort-mismatch
# hypothesis from generic Harmony oversmoothing
#
# Test A (C5): Harmony on ht70+ht71 only -- two unsorted, compositionally
#   similar samples. If the population survives strongly here, that argues
#   Harmony is fine on comparable samples, and damage in C2/C4 comes
#   specifically from forcing dissimilar (sorted) samples together.
#
# Test B (C6): Harmony on cd205neg+cd205pos only -- the two samples with the
#   most extreme, DESIGNED compositional mismatch. If the population is most
#   severely damaged here, that's direct support for the overcorrection-via-
#   sort-mismatch hypothesis.
# =============================================================================

# --- C5: ht70 + ht71 only (unsorted pair) ---
samples_unsorted_pair <- subset(samples, subset = sample.id %in% c("ht70", "ht71"))

sanity("C5 setup -- unsorted pair subset",
       "only ht70/ht71 present" = setequal(unique(samples_unsorted_pair$sample.id),
                                           c("ht70", "ht71")),
       "cells retained"         = ncol(samples_unsorted_pair) > 0
)

res5 <- run_pipeline(samples_unsorted_pair, include_ht70 = TRUE,
                     use_harmony = TRUE, label = "C5-ht70ht71only-Harmony")

# --- C6: cd205neg + cd205pos only (sorted pair, max compositional mismatch) ---
samples_sorted_pair <- subset(samples, subset = sample.id %in%
                                c("ht67-cd205neg", "ht67-cd205pos"))

sanity("C6 setup -- sorted pair subset",
       "only cd205neg/cd205pos present" = setequal(unique(samples_sorted_pair$sample.id),
                                                   c("ht67-cd205neg", "ht67-cd205pos")),
       "cells retained"                 = ncol(samples_sorted_pair) > 0
)

res6 <- run_pipeline(samples_sorted_pair, include_ht70 = TRUE,
                     use_harmony = TRUE, label = "C6-cd205pair-Harmony")

all_res <- list(res1, res2, res3, res4, res5, res6)   # <-- extend the list used in Part D


# =============================================================================
# PART C.5 -- MANUAL CLUSTER ANNOTATION (per condition, cluster numbers vary)
# =============================================================================

# Step 1: inspect one condition at a time. Start with res1.
DimPlot_scCustom(res1$object, label = TRUE) + labs(title = res1$label)
print(res1$cluster_summary)   # cTEC/mTEC-I scores + double-positive rate per cluster

# Step 2: once you've identified cluster identities by eye/markers, fill in
# the mapping below. Cluster numbers are whatever seurat_clusters shows for
# THIS condition -- will differ from other conditions.
c1_cluster_ids <- c(
  "0"  = "cTEC",
  "1"  = "cTEC",
  # ... fill in every cluster present in res1$object ...
  "8"  = "mTEC I"        # <- or whatever number is the cluster-8-equivalent
)

res1$object$cell_type_manual <- plyr::mapvalues(
  as.character(res1$object$seurat_clusters),
  from = names(c1_cluster_ids),
  to   = c1_cluster_ids
)

# Sanity check -- confirm every cluster got mapped, nothing fell through as NA
stopifnot(!any(is.na(res1$object$cell_type_manual)))
table(res1$object$cell_type_manual)

# Flag which numeric cluster you've identified as the cluster-8 equivalent
c1_c8_cluster <- "8"   # <- set this per condition once identified

DimPlot_scCustom(res1$object,
                 cells.highlight = WhichCells(res1$object, idents = c1_c8_cluster)) +
  labs(title = paste0(res1$label, ": cluster-8 equivalent = cluster ", c1_c8_cluster))

# ---- Repeat the same 3 steps for res2, res3, res4 ----




# =============================================================================
# PART D -- COMPARE CONDITIONS
# =============================================================================
cat("\n\n", strrep("=", 70), "\nPART D: CROSS-CONDITION COMPARISON\n", strrep("=", 70), "\n")

comparison <- do.call(rbind, lapply(all_res, function(r) data.frame(
  condition            = r$label,
  n_cells              = r$n_cells,
  n_clusters           = r$n_clusters,
  top_C8_cluster       = r$top_c8_cluster,
  top_C8_n_cells       = r$top_c8_n,
  top_C8_signature     = r$top_c8_sig,
  top_C8_pct_dblpos    = r$top_c8_dblpos,
  dataset_pct_dblpos   = r$overall_dblpos,
  stringsAsFactors = FALSE
)))

print(comparison)
write.csv(comparison, out("CONDITION-COMPARISON.csv"), row.names = FALSE)

cat("\n", strrep("=", 70), "\n")
cat("HOW TO INTERPRET:\n\n")
cat("A cluster-8-like population PERSISTS if, in a given condition:\n")
cat("  - some cluster still scores high on top_C8_signature, AND\n")
cat("  - that cluster's double-positive rate stays well above the dataset-wide\n")
cat("    rate (dataset_pct_dblpos) -- i.e. co-expression is still concentrated\n")
cat("    rather than smeared uniformly across all clusters.\n\n")
cat("Condition-specific readings:\n")
cat("  C1 vs C2 -> does Harmony integration dissolve it?\n")
cat("  C1 vs C3 -> is it HT70-dependent (sex / fixation-time / donor effect)?\n")
cat("  C4       -> does it survive BOTH controls simultaneously?\n\n")
cat("If the population survives C4, that is the strongest evidence available\n")
cat("from this dataset that it is not a batch or single-donor artifact.\n")
cat("If it vanishes only in C3/C4, HT70 is load-bearing -- and the 10 h vs 16 h\n")
cat("fixation difference becomes the leading candidate explanation to report.\n")
cat(strrep("=", 70), "\n")

# Side-by-side UMAPs, all four conditions
p_all <- wrap_plots(lapply(all_res, function(r) {
  DimPlot_scCustom(r$object, label = TRUE) + NoLegend() +
    labs(title = r$label) + theme(plot.title = element_text(size = 10))
}), ncol = 2)
ggsave(out("ALL-CONDITIONS-umap.pdf"), p_all, width = 14, height = 12)

# Signature localisation across conditions
p_sig <- wrap_plots(lapply(all_res, function(r) {
  FeaturePlot_scCustom(r$object, features = "C8_nearcTEC_score",
                       colors_use = viridis_light_high, order = TRUE,
                       na_cutoff = NA) +
    labs(title = r$label) + theme(plot.title = element_text(size = 10))
}), ncol = 2)
ggsave(out("ALL-CONDITIONS-C8-signature.pdf"), p_sig, width = 14, height = 12)

cat("\nAll outputs written to:", OUTDIR, "with prefix", DATE_PREFIX, "\n")