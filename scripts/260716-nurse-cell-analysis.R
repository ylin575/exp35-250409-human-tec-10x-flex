# =============================================================================
# Nurse cell subclustering -- Reviewer 2 comment on heterogeneity
#
# Goal: subcluster the nurse population, characterize subclusters, and decide
# from the data whether to claim (a) "biologically distinct subpopulations
# resolved" or (b) "heterogeneity noted but not confidently resolvable at
# current cell number/depth". Both are defensible answers -- the data decides.
#
# Key diagnostic: are subcluster top markers dominated by cell-cycle,
# ribosomal, or mitochondrial genes? -> noise-limited answer.
# Are they clean, interpretable biology (structural, engulfment, maturation)?
# -> real subpopulations answer.
# =============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(scCustomize)

tec <- readRDS("data/rds-objects/260501-after-annotation-1.rds")

sample_metadata <- data.frame(
  sample.id   = c("ht67-cd205pos", "ht67-cd205neg", "ht70", "ht71"),
  age_months  = c(3, 3, 12, 48),
  sex         = c("F", "F", "M", "F"),
  institution = c("Mount Sinai", "Mount Sinai", "Mount Sinai", "Northwell"),
  stringsAsFactors = FALSE
)

tec@meta.data <- tec@meta.data %>%
  left_join(sample_metadata, by = "sample.id")
rownames(tec@meta.data) <- colnames(tec)   # restore rownames -- dplyr strips them

# Verify
table(tec$sample.id, tec$age_months, tec$sex, tec$institution)

# -----------------------------------------------------------------------------
# 1. Subset to nurse cells and report basic stats
# -----------------------------------------------------------------------------
nurse <- subset(tec, subset = cell_type == "Nurse")

cat("=== Nurse cell subset summary ===\n")
cat("Total nurse cells:", ncol(nurse), "\n\n")

cat("Per-sample distribution:\n")
print(table(nurse$sample.id))

cat("\nnFeature_RNA distribution:\n")
print(summary(nurse$nFeature_RNA))

# Guidance for interpretation:
#   < 200 cells   -> subclustering results will be dominated by noise;
#                    the honest answer is "insufficient cell number"
#   200 - 800     -> subclustering will produce visually distinct clusters,
#                    but be very cautious about biological claims
#   > 800         -> credible subclustering; biology-driven claims defensible
#                    if markers are clean

# -----------------------------------------------------------------------------
# 2. Re-embed on just the nurse cells
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# 2. Wipe stale layers/reductions, then re-embed on just the nurse cells
# -----------------------------------------------------------------------------
DefaultAssay(nurse) <- "RNA"

# Drop stale scale.data, PCA, UMAP, neighbors, clustering carried over from
# the parent object -- keep only counts and data. Fresh pipeline runs on
# just the nurse-relevant gene set.
nurse <- DietSeurat(
  nurse,
  layers    = c("counts", "data"),
  dimreducs = NULL,
  graphs    = NULL,
  misc      = FALSE
)

nurse <- nurse %>%
  NormalizeData(verbose = TRUE) %>%
  FindVariableFeatures(nfeatures = 2000, verbose = TRUE) %>%
  ScaleData(verbose = TRUE) %>%
  RunPCA(npcs = 30, verbose = TRUE)

# Inspect how many PCs actually carry signal -- for small subsets, using too
# many PCs manufactures noise clusters
ElbowPlot(nurse, ndims = 30) + labs(title = "Nurse-only PCA: elbow plot")

# ADJUST based on elbow: for small subsets (< 1000 cells), typically 10-15 is
# enough; more risks overfitting to noise
N_PCS <- 20

nurse <- nurse %>%
  FindNeighbors(dims = 1:N_PCS, verbose = FALSE) %>%
  RunUMAP(dims = 1:N_PCS, verbose = FALSE)

# -----------------------------------------------------------------------------
# 3. Cluster at multiple resolutions -- lets you pick a defensible level
# -----------------------------------------------------------------------------
resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0)
for (r in resolutions) {
  nurse <- FindClusters(nurse, resolution = r, verbose = FALSE)
}

# What resolution produces how many subclusters?
cat("\n=== Cluster counts at each resolution ===\n")
for (r in resolutions) {
  col <- paste0("RNA_snn_res.", r)
  cat(sprintf("res %.1f  ->  %d clusters\n", r, length(unique(nurse[[col]][[1]]))))
}

# Visualize all resolutions side by side to pick one
p_res <- wrap_plots(
  lapply(resolutions, function(r) {
    col <- paste0("RNA_snn_res.", r)
    DimPlot_scCustom(nurse, group.by = col, label = TRUE, repel = TRUE) +
      labs(title = paste0("res = ", r)) +
      NoLegend()
  }),
  ncol = 3
)
p_res
ggsave("nurse_subclusters_resolutions.pdf", p_res, width = 12, height = 8)

# CHOOSE A RESOLUTION -- adjust based on the plots above.
# Rule of thumb: pick the lowest resolution that produces subclusters with
# visibly distinct positions on the UMAP.
CHOSEN_RES <- 0.2                                  # <-- adjust
Idents(nurse) <- paste0("RNA_snn_res.", CHOSEN_RES)

# -----------------------------------------------------------------------------
# 4. Sanity check: are subclusters donor-driven or biology-driven?
# -----------------------------------------------------------------------------
# If a subcluster corresponds 1:1 with a single sample.id, it's likely donor
# effects, not biology. Cross-tabulate:
cat("\n=== Subcluster vs. sample.id ===\n")
print(table(Idents(nurse), nurse$sample.id))

# Also plot side-by-side:
p_check <- (
  DimPlot_scCustom(nurse, group.by = paste0("RNA_snn_res.", CHOSEN_RES),
                   label = TRUE) + labs(title = "Nurse subclusters")
) | (
  DimPlot_scCustom(nurse, group.by = "sample.id") + labs(title = "By donor/sort")
)
p_check
ggsave("nurse_subclusters_vs_donor.pdf", p_check, width = 12, height = 5)

# -----------------------------------------------------------------------------
# 5. Find subcluster markers
# -----------------------------------------------------------------------------
nurse_markers <- FindAllMarkers(
  nurse,
  only.pos       = TRUE,
  min.pct        = 0.25,
  logfc.threshold = 0.25,
  verbose        = FALSE
)

top_markers <- nurse_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 15) %>%
  ungroup()

write.csv(nurse_markers, "nurse_subcluster_markers_ALL.csv", row.names = FALSE)
write.csv(top_markers,   "nurse_subcluster_markers_TOP15.csv", row.names = FALSE)

cat("\n=== Top markers per subcluster ===\n")
print(top_markers %>% select(cluster, gene, avg_log2FC, pct.1, pct.2, p_val_adj))

# -----------------------------------------------------------------------------
# 6. DIAGNOSTIC: is the top-marker list dominated by NOISE genes?
# -----------------------------------------------------------------------------
# Ribosomal (RPS/RPL), mitochondrial (MT-), heat-shock (HSP), and cell-cycle
# (MKI67, TOP2A, CDK1, MCM*, CCN*) genes topping the list is a strong signal
# that clustering is noise-driven rather than biology-driven.

noise_patterns <- c("^RPS", "^RPL", "^MT-", "^HSP", "^HIST",
                    "^MKI67$", "^TOP2A$", "^CDK1$", "^MCM", "^CCN")
noise_regex <- paste(noise_patterns, collapse = "|")

diagnostic <- top_markers %>%
  mutate(is_noise = grepl(noise_regex, gene)) %>%
  group_by(cluster) %>%
  summarise(
    n_top       = n(),
    n_noise     = sum(is_noise),
    pct_noise   = round(100 * n_noise / n_top, 1),
    noise_genes = paste(gene[is_noise], collapse = ", "),
    .groups = "drop"
  )

cat("\n=== Diagnostic: proportion of top markers that are noise genes ===\n")
print(diagnostic)

cat("\nINTERPRETATION:\n")
cat("  pct_noise > 40%  ->  noise-limited; frame response as 'heterogeneity\n")
cat("                       noted but not confidently resolvable at current depth'\n")
cat("  pct_noise < 20%  ->  markers look biology-driven; frame response as\n")
cat("                       'subclustering revealed distinct subpopulations\n")
cat("                       characterized by [markers]'\n")
cat("  20-40%           ->  mixed; use judgment on individual clusters\n\n")

# -----------------------------------------------------------------------------
# 7. Visualize top markers on the subcluster UMAP
# -----------------------------------------------------------------------------
# Dotplot of top markers across subclusters -- fastest way to see structure
top5 <- top_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5) %>%
  pull(gene) %>%
  unique()

p_dot <- DotPlot_scCustom(nurse, features = top5,
                          group.by = paste0("RNA_snn_res.", CHOSEN_RES)) +
  RotatedAxis() +
  labs(title = "Top-5 markers per nurse subcluster")
p_dot
ggsave("nurse_subcluster_marker_dotplot.pdf", p_dot,
       width = max(8, length(top5) * 0.35), height = 5)

# -----------------------------------------------------------------------------
# 8. Save the subclustered object
# -----------------------------------------------------------------------------
saveRDS(nurse, "nurse_subclustered.rds")

cat("\n=== Outputs ===\n")
cat("  nurse_subclusters_resolutions.pdf     -- pick a resolution\n")
cat("  nurse_subclusters_vs_donor.pdf        -- donor vs biology check\n")
cat("  nurse_subcluster_markers_ALL.csv      -- full marker table\n")
cat("  nurse_subcluster_markers_TOP15.csv    -- top 15 per subcluster\n")
cat("  nurse_subcluster_marker_dotplot.pdf   -- visual summary\n")
cat("  nurse_subclustered.rds                -- object for further analysis\n")