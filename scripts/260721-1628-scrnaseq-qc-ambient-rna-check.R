# =============================================================================
# SoupX ambient RNA check -- final artifact test for cluster 8 co-expression
#
# Ambient RNA ("soup") is estimated PER SAMPLE (each GEM well has its own
# background profile from lysed/dying cells), using the raw (unfiltered,
# includes empty droplets) + filtered (called cells only) matrices together.
#
# Logic: if PSMB11/PRSS16/LY75/KRT14/CCL19 are highly abundant in the
# ambient soup profile generally, that would suggest at least some of the
# apparent co-expression reflects background contamination rather than
# genuine per-cell signal. We test this by re-checking co-expression rates
# on SoupX-corrected counts.
# =============================================================================

library(Seurat)
library(SoupX)
library(dplyr)
library(Matrix)

DATE_PREFIX <- format(Sys.time(), "%y%m%d-%H%M")

clust8 <- readRDS("data/rds-objects/260721-1630-cluster8-artifact-checked.rds")

cTEC_markers   <- c("PSMB11", "PRSS16", "LY75")
mTEC_I_markers <- c("KRT14", "CCL19")

# -----------------------------------------------------------------------------
# 1. Paths to raw (unfiltered) and filtered .h5 files, per sample
# -----------------------------------------------------------------------------
# ADJUST to your actual raw file paths -- likely siblings of the filtered
# files you originally loaded (e.g. "sample_raw_feature_bc_matrix.h5" next
# to "sample_filtered_feature_bc_matrix.h5")

sample_paths <- list(
  "ht67-cd205neg" = list(
    raw      = "rawdata/Meyer_YL03_10x-flex-data/count/Meyer_YL03_HT67_CD205_neg/multi/count/raw_feature_bc_matrix.h5",
    filtered = "rawdata/h5-files/ht67-cd205neg/sample_filtered_feature_bc_matrix.h5"
  ),
  "ht67-cd205pos" = list(
    raw      = "rawdata/Meyer_YL03_10x-flex-data/count/Meyer_YL03_HT67_CD205_pos/multi/count/raw_feature_bc_matrix.h5",
    filtered = "rawdata/h5-files/ht67-cd205pos/sample_filtered_feature_bc_matrix.h5"
  ),
  "ht70" = list(
    raw      = "rawdata/Meyer_YL03_10x-flex-data/count/Meyer_YL03_HT70/multi/count/raw_feature_bc_matrix.h5",
    filtered = "rawdata/h5-files/ht70/sample_filtered_feature_bc_matrix.h5"
  ),
  "ht71" = list(
    raw      = "rawdata/Meyer_YL03_10x-flex-data/count/Meyer_YL03_HT71/multi/count/raw_feature_bc_matrix.h5",
    filtered = "rawdata/h5-files/ht71/sample_filtered_feature_bc_matrix.h5"
  )
)

# -----------------------------------------------------------------------------
# 2. Run SoupX per sample -- estimate + correct ambient contamination
# -----------------------------------------------------------------------------
# We only NEED corrected counts for cluster-8 cells, but SoupX's contamination
# estimation (autoEstCont) needs cluster structure across the WHOLE sample to
# find good marker genes -- so we cluster each full sample first, then extract
# corrected counts for just the cells we care about at the end.

run_soupx_for_sample <- function(sample_name, paths) {
  cat("\n=== Running SoupX for", sample_name, "===\n")
  
  raw_mat  <- Read10X_h5(paths$raw,      use.names = TRUE, unique.features = TRUE)
  filt_mat <- Read10X_h5(paths$filtered, use.names = TRUE, unique.features = TRUE)
  
  # Force identical gene set, identical order, on BOTH matrices
  common_genes <- intersect(rownames(raw_mat), rownames(filt_mat))
  cat("Raw genes:", nrow(raw_mat), " | Filtered genes:", nrow(filt_mat),
      " | Common:", length(common_genes), "\n")
  
  raw_mat  <- raw_mat[common_genes, ]
  filt_mat <- filt_mat[common_genes, ]
  
  # Quick Seurat clustering of the filtered data -- SoupX needs cluster
  # labels to identify marker genes for contamination estimation
  srat <- CreateSeuratObject(counts = filt_mat, min.cells = 3, min.features = 200)
  srat <- NormalizeData(srat, verbose = TRUE)
  srat <- FindVariableFeatures(srat, verbose = TRUE)
  srat <- ScaleData(srat, verbose = TRUE)
  srat <- RunPCA(srat, verbose = TRUE)
  srat <- FindNeighbors(srat, dims = 1:30, verbose = TRUE)
  srat <- FindClusters(srat, verbose = TRUE)
  srat <- RunUMAP(srat, dims = 1:30, verbose = TRUE)
  
  sc <- SoupChannel(raw_mat, filt_mat)
  sc <- setClusters(sc, setNames(as.character(Idents(srat)), colnames(srat)))
  sc <- setDR(sc, Embeddings(srat, "umap")[colnames(srat), ])
  
  sc <- autoEstCont(sc, verbose = TRUE)
  cat(sample_name, "estimated contamination fraction (rho):",
      round(sc$metaData$rho[1], 4), "\n")
  
  corrected_counts <- adjustCounts(sc, verbose = 1)
  
  list(
    sample_name = sample_name,
    rho = sc$metaData$rho[1],
    corrected_counts = corrected_counts,
    original_counts = filt_mat
  )
}

soupx_results <- lapply(names(sample_paths), function(s) {
  run_soupx_for_sample(s, sample_paths[[s]])
})
names(soupx_results) <- names(sample_paths)

# -----------------------------------------------------------------------------
# 3. Report contamination fraction per sample -- context for interpretation
# -----------------------------------------------------------------------------
rho_summary <- tibble(
  sample = names(soupx_results),
  contamination_fraction = sapply(soupx_results, function(x) round(x$rho, 4))
)
cat("\n=== Estimated ambient RNA contamination fraction per sample ===\n")
print(rho_summary)
write.csv(rho_summary, paste0(DATE_PREFIX, "-soupx-contamination-fraction-per-sample.csv"),
          row.names = FALSE)

cat("\nTypical contamination fractions in healthy tissue range ~1-10%. Higher\n")
cat("values (>15-20%) suggest a sample with more cell lysis/stress and\n")
cat("warrant closer scrutiny of any marginal marker-based findings.\n\n")

# -----------------------------------------------------------------------------
# 4. Extract corrected counts for cluster-8 cells specifically, re-check
#    co-expression on AMBIENT-CORRECTED data
# -----------------------------------------------------------------------------
clust8_barcodes <- colnames(clust8)   # includes sample-prefixed cell IDs

# Build a combined corrected-counts matrix across samples, matched to
# cluster 8's cell barcodes. Barcode format must match what's in clust8 --
# check and adjust the stripping logic below against your actual barcode
# naming convention (this mirrors the sample_sample-id_barcode pattern used
# elsewhere in this project).

# IMPORTANT: the original merge() call had NO add.cell.ids, so Seurat
# appended a POSITIONAL numeric suffix ("_1","_2","_3","_4") based on merge
# order, NOT a sample-name prefix. Confirmed from actual barcode formats:
#   clust8 cell names:      "AAACAAGCATCACTCG-1_1"   (bare barcode + "_N")
#   SoupX corrected_counts: "AAACAAGCAACCATTC-1"     (bare barcode, no suffix)
# Merge order from the original script:
#   merge(x = ht67_cd205neg_sub, y = list(ht67_cd205pos_sub, ht70_sub, ht71_sub))
#   -> suffix _1 = ht67-cd205neg, _2 = ht67-cd205pos, _3 = ht70, _4 = ht71

merge_order_suffix <- c(
  "ht67-cd205neg" = "_1",
  "ht67-cd205pos" = "_2",
  "ht70"          = "_3",
  "ht71"          = "_4"
)

get_corrected_for_barcodes <- function(sample_name, result, target_barcodes) {
  suffix <- merge_order_suffix[[sample_name]]
  # Keep only this sample's cells, then strip the positional suffix to get
  # the bare barcode that matches SoupX's raw per-sample matrix column names
  this_sample_barcodes <- target_barcodes[endsWith(target_barcodes, suffix)]
  bare <- sub(paste0(suffix, "$"), "", this_sample_barcodes)
  bare <- intersect(bare, colnames(result$corrected_counts))
  if (length(bare) == 0) return(NULL)
  mat <- result$corrected_counts[, bare, drop = FALSE]
  # Re-attach the original clust8-style names (with suffix) as column names
  # so this matrix can be cbind'd cleanly and matched back later if needed
  colnames(mat) <- paste0(bare, suffix)
  mat
}

corrected_list <- lapply(names(soupx_results), function(s) {
  get_corrected_for_barcodes(s, soupx_results[[s]], clust8_barcodes)
})
corrected_list <- corrected_list[!sapply(corrected_list, is.null)]

cat("Matched cells per sample:\n")
print(sapply(corrected_list, ncol))

# Combine into one matrix (genes as rows, cluster-8 cells as columns)
common_genes <- Reduce(intersect, lapply(corrected_list, rownames))
corrected_combined <- do.call(cbind, lapply(corrected_list, function(m) m[common_genes, ]))

cat("Corrected counts matched for", ncol(corrected_combined),
    "of", ncol(clust8), "cluster 8 cells.\n")
cat("(If this is much lower than 1300, check the barcode-matching logic above\n")
cat("against your actual naming convention before trusting the result.)\n\n")

# -----------------------------------------------------------------------------
# 5. Re-test co-expression on ambient-CORRECTED raw counts
# -----------------------------------------------------------------------------
MIN_UMI <- 2

cTEC_markers_avail   <- intersect(cTEC_markers, rownames(corrected_combined))
mTEC_I_markers_avail <- intersect(mTEC_I_markers, rownames(corrected_combined))

cTEC_corrected_sum   <- Matrix::colSums(corrected_combined[cTEC_markers_avail, , drop = FALSE])
mTEC_I_corrected_sum <- Matrix::colSums(corrected_combined[mTEC_I_markers_avail, , drop = FALSE])

corrected_coexp <- tibble(
  cell = colnames(corrected_combined),
  cTEC_corrected_pos   = cTEC_corrected_sum >= MIN_UMI,
  mTEC_I_corrected_pos = mTEC_I_corrected_sum >= MIN_UMI
) %>%
  dplyr::count(cTEC_corrected_pos, mTEC_I_corrected_pos) %>%
  mutate(pct = round(100 * n / sum(n), 1))

cat("=== Co-expression rate on AMBIENT-CORRECTED counts (>=", MIN_UMI, "UMI) ===\n")
print(corrected_coexp)

write.csv(corrected_coexp, paste0(DATE_PREFIX, "-cluster8-ambient-corrected-coexpression.csv"),
          row.names = FALSE)

cat("\nCOMPARE this double-positive percentage to the pre-correction strict\n")
cat("threshold result (87.4%). If the corrected rate is similar, ambient RNA\n")
cat("is NOT a meaningful driver of the co-expression finding. A substantial\n")
cat("drop after correction would indicate some of the signal was ambient\n")
cat("background rather than genuine per-cell expression.\n")

# -----------------------------------------------------------------------------
# 6. Save
# -----------------------------------------------------------------------------
saveRDS(soupx_results, paste0(DATE_PREFIX, "-soupx-results-all-samples.rds"))
