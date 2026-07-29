# =============================================================================
# Ragazzini et al. 2023 PolyKRT/TESC module score
#
# Full published marker set (Dev Cell 2023, Figure 1B legend), CTGF excluded
# per user (not found in probe panel; checked for renamed alias CCN2 as
# fallback before concluding true absence).
#
# Tests Reviewer 2's specific request: does this externally-published
# TESC/immature-TEC signature concentrate into a distinct cluster in OUR
# data, or spread across subsets?
# =============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(scCustomize)

DATE_PREFIX <- format(Sys.time(), "%y%m%d-%H%M")

tec <- readRDS("data/rds-objects/260501-after-annotation-1.rds")

# -----------------------------------------------------------------------------
# 1. Define the module -- full Ragazzini PolyKRT signature, minus CTGF
# -----------------------------------------------------------------------------
polyKRT_module_genes <- c(
  "KRT13", "KRT14", "KRT15", "KRT17", "KRT19",
  "CCL19", "CEBPD", "CLU", "FN1", "IFITM3",
  "TIMP1", "VCAM1", "TAGLN", "BCAM", "LIFR", "CH25H", "CCN2"
)

# Check for CTGF under its current HGNC symbol (CCN2) before concluding
# it's genuinely absent from the probe panel
if ("CCN2" %in% rownames(tec) && !("CTGF" %in% rownames(tec))) {
  cat("NOTE: CTGF not found, but CCN2 (its current HGNC symbol) IS present.\n")
  cat("      Not auto-added (per instruction to exclude CTGF) -- add 'CCN2'\n")
  cat("      to polyKRT_module_genes above if you want it included.\n\n")
} else if (!("CTGF" %in% rownames(tec)) && !("CCN2" %in% rownames(tec))) {
  cat("NOTE: Neither CTGF nor renamed symbol CCN2 found. Confirmed absent.\n\n")
}

missing <- setdiff(polyKRT_module_genes, rownames(tec))
if (length(missing) > 0) {
  warning("Not in object, dropping from module: ", paste(missing, collapse = ", "))
  polyKRT_module_genes <- setdiff(polyKRT_module_genes, missing)
}
cat("Module genes used (n=", length(polyKRT_module_genes), "):\n", sep = "")
cat(paste(polyKRT_module_genes, collapse = ", "), "\n\n")

# -----------------------------------------------------------------------------
# 2. Run AddModuleScore
# -----------------------------------------------------------------------------
tec <- AddModuleScore(
  tec,
  features = list(polyKRT_module_genes),
  name     = "PolyKRT_module",
  verbose  = TRUE
)

tec$PolyKRT_module_score <- tec$PolyKRT_module1
tec$PolyKRT_module1 <- NULL

cat("Score range:", round(range(tec$PolyKRT_module_score), 3), "\n\n")

# -----------------------------------------------------------------------------
# 3. Visualize
# -----------------------------------------------------------------------------
p_vln <- VlnPlot_scCustom(
  tec, features = "PolyKRT_module_score", group.by = "cell_type", pt.size = 0.1
) +
  labs(title = "Ragazzini PolyKRT/TESC module score by cell type",
       y = "Module score")
print(p_vln)
ggsave(paste0(DATE_PREFIX, "-PolyKRT-module-score-violin.pdf"), p_vln,
       width = 8, height = 5)

p_umap <- FeaturePlot_scCustom(
  tec, features = "PolyKRT_module_score", colors_use = viridis_light_high,
  order = TRUE, pt.size = 0.4
) +
  labs(title = "PolyKRT/TESC module score, UMAP")
print(p_umap)
ggsave(paste0(DATE_PREFIX, "-PolyKRT-module-score-umap.pdf"), p_umap,
       width = 6, height = 5)

# -----------------------------------------------------------------------------
# 4. QUANTITATIVE TEST: mean score per cell type + pairwise significance
# -----------------------------------------------------------------------------
cat("=== Mean PolyKRT module score per cell type ===\n")
score_summary <- tec@meta.data %>%
  group_by(cell_type) %>%
  summarise(
    n            = n(),
    mean_score   = round(mean(PolyKRT_module_score), 3),
    median_score = round(median(PolyKRT_module_score), 3),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_score))
print(score_summary)

write.csv(score_summary, paste0(DATE_PREFIX, "-PolyKRT-module-score-by-celltype.csv"),
          row.names = FALSE)

top_type <- score_summary$cell_type[1]
other_types <- setdiff(unique(tec$cell_type), top_type)

pairwise_tests <- lapply(other_types, function(ct) {
  top_scores   <- tec$PolyKRT_module_score[tec$cell_type == top_type]
  other_scores <- tec$PolyKRT_module_score[tec$cell_type == ct]
  test <- wilcox.test(top_scores, other_scores)
  tibble(
    comparison = paste0(top_type, " vs ", ct),
    p_value    = test$p.value,
    top_mean   = mean(top_scores),
    other_mean = mean(other_scores)
  )
}) %>% bind_rows() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

cat("\n=== Highest-scoring cell type (", top_type,
    ") vs. every other type (Wilcoxon, BH-adjusted) ===\n", sep = "")
print(pairwise_tests)

write.csv(pairwise_tests, paste0(DATE_PREFIX, "-PolyKRT-module-pairwise-tests.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 5. CONCENTRATION TEST
# -----------------------------------------------------------------------------
N_TOP <- score_summary$n[1]

top_scoring <- tec@meta.data %>%
  slice_max(order_by = PolyKRT_module_score, n = N_TOP)

concentration <- top_scoring %>%
  count(cell_type, sort = TRUE) %>%
  mutate(pct_of_top_group = round(100 * n / N_TOP, 1))

cat("\n=== Cell type composition of the top", N_TOP,
    "highest PolyKRT-scoring cells ===\n")
print(concentration)

write.csv(concentration, paste0(DATE_PREFIX, "-PolyKRT-module-concentration-test.csv"),
          row.names = FALSE)

cat("\nINTERPRETATION:\n")
cat("  If top scorers concentrate heavily (e.g. >70%) in ONE existing cell\n")
cat("  type, the Ragazzini PolyKRT/TESC signature DOES map to a distinct,\n")
cat("  identifiable population in our data.\n")
cat("  If top scorers spread broadly across MULTIPLE cell types with no\n")
cat("  single dominant group, this supports the position that this\n")
cat("  signature does not cleanly demarcate a distinct TESC/immature-TEC\n")
cat("  population in our dataset -- report honestly either way.\n\n")

# -----------------------------------------------------------------------------
# 6. Does PolyKRT signature overlap with BMC-TEC specifically?
# -----------------------------------------------------------------------------
if ("BMC_module_score" %in% colnames(tec@meta.data)) {
  cat("=== Correlation: PolyKRT score vs. your own BMC-TEC module score ===\n")
  corr <- cor.test(tec$PolyKRT_module_score, tec$BMC_module_score,
                   method = "spearman")
  cat("Spearman rho =", round(corr$estimate, 3),
      ", p =", format.pval(corr$p.value), "\n\n")
  
  p_corr <- ggplot(tec@meta.data,
                   aes(x = BMC_module_score, y = PolyKRT_module_score,
                       color = cell_type)) +
    geom_point(size = 0.5, alpha = 0.5) +
    labs(title = "BMC-TEC module vs. PolyKRT/TESC module, per cell",
         x = "BMC-TEC module score", y = "PolyKRT/TESC module score") +
    theme_minimal()
  print(p_corr)
  ggsave(paste0(DATE_PREFIX, "-BMC-vs-PolyKRT-module-correlation.pdf"), p_corr,
         width = 7, height = 5)
} else {
  cat("NOTE: BMC_module_score not found -- run the BMC-TEC module script\n")
  cat("      first if you want the cross-module comparison in Section 6.\n")
}

# -----------------------------------------------------------------------------
# 7. Dotplot of individual module genes + composite score, by cell type
# -----------------------------------------------------------------------------
p_dot <- DotPlot_scCustom(
  tec, features = c(polyKRT_module_genes, "PolyKRT_module_score"),
  group.by = "cell_type"
) +
  RotatedAxis() +
  labs(title = "PolyKRT module genes + composite score, by cell type")
print(p_dot)
ggsave(paste0(DATE_PREFIX, "-PolyKRT-module-genes-dotplot.pdf"), p_dot,
       width = 10, height = 5)

# -----------------------------------------------------------------------------
# 8. Save
# -----------------------------------------------------------------------------
saveRDS(tec, paste0(DATE_PREFIX, "-tec-with-PolyKRT-module-score.rds"))
