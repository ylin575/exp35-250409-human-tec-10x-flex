# =============================================================================
# CD49f/BCAM module score -- Reviewer 2's request
#
# "the clusters the authors were following with CD49f and BCAM are also
# defined at the RNA level by multiple genes that the authors should gather
# into a module and score across their clusters, to test whether these
# populations form a distinct cluster or are spread across subsets."
# =============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(scCustomize)

DATE_PREFIX <- format(Sys.time(), "%y%m%d-%H%M")

tec <- readRDS("data/rds-objects/260501-after-annotation-1.rds")

# -----------------------------------------------------------------------------
# 1. Define the module -- BMC-TEC identity genes (CD49f/BCAM + ECM block)
# -----------------------------------------------------------------------------
# CD200 deliberately excluded: BMC-TEC is CD200-LOW, so it's a negative
# marker for this population, not a positive one -- including it in an
# additive score would dilute rather than sharpen the signal.

bmc_module_genes <- c("COL4A5", "COL4A6", "ITGB4", "COL8A1")

missing <- setdiff(bmc_module_genes, rownames(tec))
if (length(missing) > 0) {
  warning("Not in object, dropping from module: ", paste(missing, collapse = ", "))
  bmc_module_genes <- setdiff(bmc_module_genes, missing)
}
cat("Module genes used:", paste(bmc_module_genes, collapse = ", "), "\n\n")

# -----------------------------------------------------------------------------
# 2. Run AddModuleScore
# -----------------------------------------------------------------------------
# features MUST be a list (even for one module) -- AddModuleScore is built to
# score multiple modules in one call. name= gets a number appended (e.g.
# "BMC_module1"), NOT "BMC_module" -- rename immediately to avoid confusion
# in every downstream line.

tec <- AddModuleScore(
  tec,
  features = list(bmc_module_genes),
  name     = "BMC_module",
  verbose  = TRUE
)

# Fix the auto-appended "1" suffix
tec$BMC_module_score <- tec$BMC_module1
tec$BMC_module1 <- NULL

cat("Score range:", round(range(tec$BMC_module_score), 3), "\n\n")

# -----------------------------------------------------------------------------
# 3. Visualize -- violin plot across all cell types (the direct visual answer)
# -----------------------------------------------------------------------------
p_vln <- VlnPlot_scCustom(
  tec, features = "BMC_module_score", group.by = "cell_type", pt.size = 0.1
) +
  labs(title = "CD49f/BCAM module score by cell type",
       y = "Module score")
print(p_vln)
ggsave(paste0(DATE_PREFIX, "-BMC-module-score-violin.pdf"), p_vln,
       width = 8, height = 5)

# UMAP feature plot -- shows spatial concentration (or lack thereof)
p_umap <- FeaturePlot_scCustom(
  tec, features = "BMC_module_score", colors_use = viridis_light_high,
  order = TRUE, pt.size = 0.4
) +
  labs(title = "CD49f/BCAM module score, UMAP")
print(p_umap)
ggsave(paste0(DATE_PREFIX, "-BMC-module-score-umap.pdf"), p_umap,
       width = 6, height = 5)

# -----------------------------------------------------------------------------
# 4. QUANTITATIVE TEST: does BMC-TEC score significantly higher than every
#    other cell type? (directly answers "distinct cluster" vs "spread")
# -----------------------------------------------------------------------------
cat("=== Mean module score per cell type ===\n")
score_summary <- tec@meta.data %>%
  group_by(cell_type) %>%
  summarise(
    n         = n(),
    mean_score = round(mean(BMC_module_score), 3),
    median_score = round(median(BMC_module_score), 3),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_score))
print(score_summary)

write.csv(score_summary, paste0(DATE_PREFIX, "-BMC-module-score-by-celltype.csv"),
          row.names = FALSE)

# Pairwise Wilcoxon: BMC-TEC vs every other cell type, individually
other_types <- setdiff(unique(tec$cell_type), "BMC+ TEC")

pairwise_tests <- lapply(other_types, function(ct) {
  bmc_scores   <- tec$BMC_module_score[tec$cell_type == "BMC+ TEC"]
  other_scores <- tec$BMC_module_score[tec$cell_type == ct]
  test <- wilcox.test(bmc_scores, other_scores)
  tibble(
    comparison = paste0("BMC+ TEC vs ", ct),
    p_value    = test$p.value,
    bmc_mean   = mean(bmc_scores),
    other_mean = mean(other_scores)
  )
}) %>% bind_rows() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

cat("\n=== BMC+ TEC vs. every other cell type (Wilcoxon, BH-adjusted) ===\n")
print(pairwise_tests)

write.csv(pairwise_tests, paste0(DATE_PREFIX, "-BMC-module-pairwise-tests.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 5. CONCENTRATION TEST: of the top-scoring cells dataset-wide, what fraction
#    actually falls within the BMC-TEC cluster vs. elsewhere?
#    This directly operationalizes "distinct cluster" vs "spread across
#    subsets" as a single number.
# -----------------------------------------------------------------------------
n_bmc <- sum(tec$cell_type == "BMC+ TEC")

top_scoring <- tec@meta.data %>%
  slice_max(order_by = BMC_module_score, n = n_bmc)   # same N as BMC-TEC cluster size

concentration <- top_scoring %>%
  count(cell_type, sort = TRUE) %>%
  mutate(pct_of_top_group = round(100 * n / n_bmc, 1))

cat("\n=== Cell type composition of the top", n_bmc,
    "highest-scoring cells (matched to BMC-TEC cluster size) ===\n")
print(concentration)

write.csv(concentration, paste0(DATE_PREFIX, "-BMC-module-concentration-test.csv"),
          row.names = FALSE)

cat("\nINTERPRETATION:\n")
cat("  If the top", n_bmc, "highest-scoring cells are overwhelmingly (e.g. >80%)\n")
cat("  annotated as BMC+ TEC, the module score IS concentrated in a distinct\n")
cat("  cluster -- supports BMC-TEC as a real, RNA-defined population.\n")
cat("  If high scorers are spread substantially across other cell types too,\n")
cat("  the module does NOT cleanly define a distinct population at the RNA\n")
cat("  level -- report this honestly, it directly answers the reviewer's\n")
cat("  question either way.\n")

# -----------------------------------------------------------------------------
# 6. Dotplot showing individual module genes alongside the composite score
# -----------------------------------------------------------------------------
p_dot <- DotPlot_scCustom(
  tec, features = c(bmc_module_genes, "BMC_module_score"), group.by = "cell_type"
) +
  RotatedAxis() +
  labs(title = "Module genes + composite score, by cell type")
print(p_dot)
ggsave(paste0(DATE_PREFIX, "-BMC-module-genes-and-score-dotplot.pdf"), p_dot,
       width = 8, height = 5)

# -----------------------------------------------------------------------------
# 7. Save
# -----------------------------------------------------------------------------
saveRDS(tec, paste0(DATE_PREFIX, "-tec-with-BMC-module-score.rds"))
