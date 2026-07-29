# =============================================================================
# Standard scRNAseq artifact checks -- cluster 8 cTEC/mTEC-I co-expression
#
# This population isn't described in the literature, so before reporting it
# as real biology, systematically rule out the standard technical
# explanations for "cell expressing markers of two different cell types":
#
#   1. DOUBLETS        -- the single most likely technical explanation for
#                          apparent co-expression of two lineage markers
#   2. QC OUTLIERS      -- do double-positive cells have anomalous
#                          nCount/nFeature/percent.mt (doublets often have
#                          inflated nCount/nFeature from two cells' RNA)
#   3. RAW UMI FLOOR    -- module score >0 can be driven by very few raw
#                          reads; re-test co-expression with a stricter,
#                          raw-count-based threshold
#   4. EXPECTED RATE    -- compare observed "co-expressing" rate against
#                          the expected technical doublet rate for this
#                          platform/loading -- if co-expression rate is far
#                          higher than expected doublets, doublets alone
#                          cannot explain it
#   5. CELL CYCLE       -- rule out a cycling-state confound
#   6. AMBIENT RNA       -- optional, flagged as requiring raw unfiltered
#                          matrix (SoupX/decontX), noted but not run here
#                          unless that data is available
# =============================================================================

library(Seurat)
library(SingleCellExperiment)
library(scDblFinder)
library(dplyr)
library(ggplot2)
library(scCustomize)
library(patchwork)

DATE_PREFIX <- format(Sys.time(), "%y%m%d-%H%M")

clust8 <- readRDS("data/rds-objects/260721-1508-cluster8-subclustered.rds")

# -----------------------------------------------------------------------------
# 1. DOUBLET DETECTION -- scDblFinder, respecting sample/batch structure
# -----------------------------------------------------------------------------
# Verified current best practice (independently benchmarked top performer,
# Xi & Li 2021; scran itself now recommends scDblFinder over its own older
# doubletCells method). Doublets only form WITHIN a GEM well, so this MUST
# be run with samples= set to your original sample.id -- running across
# pooled/merged samples without this would incorrectly compare cells from
# different wells as potential doublet partners.

sce <- as.SingleCellExperiment(clust8, assay = "RNA")

set.seed(42)   # scDblFinder involves randomness (simulated doublets) -- fix for reproducibility
sce <- scDblFinder(sce, samples = "sample.id", verbose = TRUE)

clust8$scDblFinder_score <- sce$scDblFinder.score
clust8$scDblFinder_class <- sce$scDblFinder.class   # "singlet" or "doublet"

cat("=== scDblFinder classification ===\n")
print(table(clust8$scDblFinder_class))
print(round(prop.table(table(clust8$scDblFinder_class)) * 100, 1))

# -----------------------------------------------------------------------------
# 2. EXPECTED MULTIPLET RATE BENCHMARK
# -----------------------------------------------------------------------------
# 10x's published multiplet-rate scales with cells loaded/recovered.
# Standard chemistry: ~0.8% per 1000 cells recovered.
# Newer GEM-X/HT chemistry (closer to Flex): roughly half that, ~0.4%/1000.
# NOTE: exact Flex-specific rate not verified here -- check your 10x Flex
# user guide for the precise expected multiplet rate at your loading target;
# treat the estimate below as an approximate benchmark, not a precise figure.

cells_per_sample <- table(clust8$sample.id)  # this is cluster-8-only count;
# for an accurate expected rate you'd want the ORIGINAL per-sample total cell
# count before subsetting to cluster 8 -- substitute those numbers if available
cat("\n=== Approximate cells per sample (cluster 8 subset only -- see note above) ===\n")
print(cells_per_sample)

approx_expected_multiplet_pct <- 0.4 * (cells_per_sample / 1000)  # GEM-X/HT-style estimate
cat("\nApproximate expected multiplet rate per sample (rough benchmark):\n")
print(round(approx_expected_multiplet_pct, 2))

cat("\nCompare this to the co-expression rate found earlier (67.5% of cluster 8).\n")
cat("If expected doublet rate is in the low single digits (%), doublets alone\n")
cat("cannot explain a 67.5% double-positive rate -- some other explanation\n")
cat("(real biology, ambient RNA, or a QC/threshold artifact) must dominate.\n\n")

# -----------------------------------------------------------------------------
# 3. Does doublet status explain the double-positive cells specifically?
# -----------------------------------------------------------------------------
# The key cross-tab: of cells called double-positive (cTEC+ AND mTEC-I+),
# what fraction are ALSO flagged as doublets by scDblFinder?

cat("=== Doublet status vs. co-expression status ===\n")
doublet_vs_coexp <- clust8@meta.data %>%
  dplyr::count(scDblFinder_class, cTEC_pos, mTEC_I_pos) %>%
  dplyr::group_by(cTEC_pos, mTEC_I_pos) %>%
  dplyr::mutate(pct_within_group = round(100 * n / sum(n), 1)) %>%
  dplyr::ungroup()
print(doublet_vs_coexp)

write.csv(doublet_vs_coexp, paste0(DATE_PREFIX, "-cluster8-doublet-vs-coexpression.csv"),
          row.names = FALSE)

cat("\nINTERPRETATION: if the double-positive group has a MUCH higher doublet\n")
cat("rate than the single-positive or negative groups, doublets are a\n")
cat("meaningful contributor. If the doublet rate is similarly low across ALL\n")
cat("groups (consistent with the low expected background rate from Section 2),\n")
cat("doublets are NOT the primary explanation for co-expression.\n\n")

# -----------------------------------------------------------------------------
# 4. QC OUTLIER CHECK -- do double-positive cells look like doublets by
#    nCount/nFeature (independent of the formal doublet score)?
# -----------------------------------------------------------------------------
clust8$coexp_group <- case_when(
  clust8$cTEC_pos & clust8$mTEC_I_pos  ~ "Double-positive",
  clust8$cTEC_pos & !clust8$mTEC_I_pos ~ "cTEC-only",
  !clust8$cTEC_pos & clust8$mTEC_I_pos ~ "mTEC-I-only",
  TRUE ~ "Double-negative"
)

cat("=== QC metrics by co-expression group ===\n")
qc_by_group <- clust8@meta.data %>%
  group_by(coexp_group) %>%
  summarise(
    n = n(),
    mean_nCount   = round(mean(nCount_RNA), 0),
    mean_nFeature = round(mean(nFeature_RNA), 0),
    mean_pct_mt   = round(mean(percent.mt), 2),
    .groups = "drop"
  )
print(qc_by_group)

write.csv(qc_by_group, paste0(DATE_PREFIX, "-cluster8-QC-by-coexpression-group.csv"),
          row.names = FALSE)

p_qc <- VlnPlot_scCustom(clust8, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
                         group.by = "coexp_group", pt.size = 0.1, num_columns = 3)
print(p_qc)
ggsave(paste0(DATE_PREFIX, "-cluster8-QC-by-coexpression-violin.pdf"), p_qc,
       width = 12, height = 4)

cat("\nINTERPRETATION: if 'Double-positive' cells have notably HIGHER\n")
cat("nCount/nFeature than other groups, that's consistent with doublets\n")
cat("(two cells' worth of RNA). Similar nCount/nFeature across groups argues\n")
cat("against a simple doublet explanation.\n\n")

# -----------------------------------------------------------------------------
# 5. RAW UMI FLOOR -- stricter co-expression criterion using actual counts,
#    not just module score > 0 (which can be driven by very sparse signal)
# -----------------------------------------------------------------------------
cTEC_markers   <- c("PSMB11", "PRSS16", "LY75")
mTEC_I_markers <- c("KRT14", "CCL19")
cTEC_markers   <- intersect(cTEC_markers, rownames(clust8))
mTEC_I_markers <- intersect(mTEC_I_markers, rownames(clust8))

counts_mat <- GetAssayData(clust8, assay = "RNA", layer = "counts")

# Stricter rule: require at least MIN_UMI raw counts summed across the
# marker panel, not just a positive (but possibly noise-driven) module score
MIN_UMI <- 2

cTEC_raw_sum   <- Matrix::colSums(counts_mat[cTEC_markers, , drop = FALSE])
mTEC_I_raw_sum <- Matrix::colSums(counts_mat[mTEC_I_markers, , drop = FALSE])

clust8$cTEC_raw_pos   <- cTEC_raw_sum >= MIN_UMI
clust8$mTEC_I_raw_pos <- mTEC_I_raw_sum >= MIN_UMI

strict_coexp <- clust8@meta.data %>%
  dplyr::count(cTEC_raw_pos, mTEC_I_raw_pos) %>%
  mutate(pct = round(100 * n / sum(n), 1))

cat("=== STRICT co-expression (>=", MIN_UMI, "raw UMI per marker panel) ===\n")
print(strict_coexp)

strict_double_pos_pct <- strict_coexp %>%
  filter(cTEC_raw_pos, mTEC_I_raw_pos) %>%
  pull(pct)

cat("\nStrict double-positive rate:", strict_double_pos_pct, "%\n")
cat("(Compare to the module-score-based 67.5% from before -- if this stricter,\n")
cat(" raw-count-based rate is still substantial, that further argues against\n")
cat(" a threshold/noise artifact explanation.)\n\n")

write.csv(strict_coexp, paste0(DATE_PREFIX, "-cluster8-strict-raw-UMI-coexpression.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 6. CELL CYCLE CONFOUND CHECK
# -----------------------------------------------------------------------------
s_genes <- cc.genes.updated.2019$s.genes
g2m_genes <- cc.genes.updated.2019$g2m.genes
clust8 <- CellCycleScoring(clust8, s.features = s_genes, g2m.features = g2m_genes,
                           set.ident = FALSE, verbose = TRUE)

cat("=== Cell cycle phase by co-expression group ===\n")
print(table(clust8$Phase, clust8$coexp_group))

# -----------------------------------------------------------------------------
# 7. AMBIENT RNA -- flagged, not run (requires raw UNFILTERED matrix)
# -----------------------------------------------------------------------------
cat("\n", strrep("=", 70), "\n")
cat("NOT RUN: Ambient RNA contamination check (SoupX or decontX).\n")
cat("These tools require the RAW, unfiltered barcode matrix (including\n")
cat("empty droplets) to estimate the ambient RNA 'soup' profile -- not just\n")
cat("the filtered/QC'd object used here. If you have access to the raw\n")
cat("(unfiltered) Cell Ranger output per sample, this is worth running as a\n")
cat("final check: if PSMB11/PRSS16/LY75/KRT14/CCL19 are highly abundant in\n")
cat("the ambient soup profile generally (common for highly-expressed genes\n")
cat("in dying/lysed cells), that would support an ambient contamination\n")
cat("explanation rather than genuine per-cell co-expression.\n")
cat(strrep("=", 70), "\n\n")

# -----------------------------------------------------------------------------
# 8. SUMMARY
# -----------------------------------------------------------------------------
cat("=== SUMMARY: artifact checks completed ===\n")
cat("1. Doublet detection (scDblFinder):        see Section 1 + 3\n")
cat("2. Expected multiplet rate benchmark:      see Section 2\n")
cat("3. QC outlier check (nCount/nFeature):     see Section 4\n")
cat("4. Raw UMI floor (stricter co-expression): see Section 5\n")
cat("5. Cell cycle confound:                    see Section 6\n")
cat("6. Ambient RNA (requires raw matrix):      NOT RUN -- see Section 7 note\n")

saveRDS(clust8, paste0(DATE_PREFIX, "-cluster8-artifact-checked.rds"))
