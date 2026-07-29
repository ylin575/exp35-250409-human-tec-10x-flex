# =============================================================================
# Cluster 8 (mTEC I) investigation -- Reviewer 2's comment on mTEC I cells
# appearing at the bottom of the cTEC cluster, in both CD205+ and all-TEC data
#
# Observation: cluster 8, subsetted alone, splits into two spatially distinct
# UMAP regions -- one adjacent to cTEC/BMC-TEC, one in the proper mTEC I
# neighborhood. This script:
#   1. Subclusters cluster 8 objectively (not by manual UMAP-coordinate gating)
#   2. Tests cTEC vs mTEC I marker co-expression per resulting sub-population
#   3. Checks whether the cTEC-adjacent sub-population is CD205+-sort-enriched
# =============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(scCustomize)
library(patchwork)

DATE_PREFIX <- format(Sys.time(), "%y%m%d-%H%M")

tec <- readRDS("data/rds-objects/260501-after-annotation-1.rds") # the "test" 
# object, with seurat_clusters holding the original 0-26 cluster labels AND 
# cell_type holding the mTEC I / cTEC / etc annotation

# -----------------------------------------------------------------------------
# 1. Subset cluster 8 specifically (using the ORIGINAL numeric cluster ID,
#    not the collapsed cell_type annotation -- other clusters also map to
#    "mTEC I" and would dilute this analysis)
# -----------------------------------------------------------------------------
# Confirm the original per-cell cluster ID column name -- likely
# "seurat_clusters" or similar from your original script's FindClusters() call
stopifnot("seurat_clusters" %in% colnames(tec@meta.data))

clust8 <- subset(tec, subset = seurat_clusters == "8")
cat("Cluster 8 total cells:", ncol(clust8), "\n")
cat("Per-sample distribution:\n")
print(table(clust8$sample.id))

# -----------------------------------------------------------------------------
# 2. Wipe stale layers/reductions from parent, re-embed on cluster 8 alone
# -----------------------------------------------------------------------------
DefaultAssay(clust8) <- "RNA"
clust8 <- DietSeurat(clust8, layers = c("counts", "data"),
                     dimreducs = NULL, graphs = NULL, misc = FALSE)

clust8 <- clust8 %>%
  NormalizeData(verbose = TRUE) %>%
  FindVariableFeatures(nfeatures = 2000, verbose = TRUE) %>%
  ScaleData(verbose = TRUE) %>%
  RunPCA(npcs = 30, verbose = TRUE)

ElbowPlot(clust8, ndims = 30) + labs(title = "Cluster 8-only PCA: elbow plot")
ggsave(paste0(DATE_PREFIX, "-cluster8-elbow.pdf"), width = 6, height = 5)

N_PCS <- 10   # adjust based on elbow plot and cluster 8 cell count

clust8 <- clust8 %>%
  FindNeighbors(dims = 1:N_PCS, verbose = TRUE) %>%
  RunUMAP(dims = 1:N_PCS, verbose = TRUE)

# -----------------------------------------------------------------------------
# 3. Subcluster at a coarse resolution first -- looking for the 2-group split
# -----------------------------------------------------------------------------
clust8 <- FindClusters(clust8, resolution = 0.2, verbose = TRUE)

p_sub <- DimPlot_scCustom(clust8, group.by = "RNA_snn_res.0.2", label = TRUE) +
  labs(title = "Cluster 8 sub-clustering (res=0.2)")
print(p_sub)
ggsave(paste0(DATE_PREFIX, "-cluster8-subclusters.pdf"), p_sub, width = 6, height = 5)

# Sanity check: does this new sub-clustering correspond to the two spatial
# regions seen on the ORIGINAL full-dataset UMAP coordinates?
# (Requires the original umap embedding's coordinates to still be attached,
# or re-plot using the parent object's UMAP colored by these sub-cluster IDs)
clust8$subcluster <- clust8$RNA_snn_res.0.2

cat("\n=== Sub-cluster sizes ===\n")
print(table(clust8$subcluster))

# -----------------------------------------------------------------------------
# 4. Sample composition per sub-cluster -- tests the CD205+ enrichment question
# -----------------------------------------------------------------------------
cat("\n=== Sub-cluster vs. sample.id ===\n")
composition <- table(clust8$subcluster, clust8$sample.id)
print(composition)

comp_pct <- prop.table(composition, margin = 1) * 100
cat("\n=== As % within each sub-cluster ===\n")
print(round(comp_pct, 1))

write.csv(as.data.frame.matrix(composition),
          paste0(DATE_PREFIX, "-cluster8-subcluster-vs-sample.csv"))

# INTERPRETATION: if the cTEC-adjacent sub-population is disproportionately
# ht67-cd205pos relative to its baseline share, that supports a CD205-sort-
# specific explanation. If it's proportionally represented across ALL
# samples including unsorted ht70/ht71, that supports a general biological
# (co-expression/transitional) explanation independent of sorting.

# -----------------------------------------------------------------------------
# 5. cTEC vs mTEC I marker co-expression -- the direct answer to the reviewer
# -----------------------------------------------------------------------------
# ADJUST these panels if your manuscript has established a more specific
# cTEC / mTEC I marker set beyond these canonical choices.
cTEC_markers   <- c("PSMB11", "PRSS16", "LY75")     # canonical cTEC identity
mTEC_I_markers <- c("KRT14", "CCL19")               # pre-AIRE / mTEC I markers

missing_cTEC <- setdiff(cTEC_markers, rownames(clust8))
missing_mTEC <- setdiff(mTEC_I_markers, rownames(clust8))
if (length(missing_cTEC) > 0) warning("Missing cTEC markers: ", paste(missing_cTEC, collapse=", "))
if (length(missing_mTEC) > 0) warning("Missing mTEC I markers: ", paste(missing_mTEC, collapse=", "))
cTEC_markers   <- intersect(cTEC_markers, rownames(clust8))
mTEC_I_markers <- intersect(mTEC_I_markers, rownames(clust8))

clust8 <- AddModuleScore(clust8, features = list(cTEC_markers),
                         name = "cTEC_score", verbose = TRUE)
clust8 <- AddModuleScore(clust8, features = list(mTEC_I_markers),
                         name = "mTEC_I_score", verbose = TRUE)
clust8$cTEC_score   <- clust8$cTEC_score1;   clust8$cTEC_score1   <- NULL
clust8$mTEC_I_score <- clust8$mTEC_I_score1; clust8$mTEC_I_score1 <- NULL

# Scatter: does each sub-cluster occupy a distinct region of cTEC-vs-mTEC-I
# score space, or is there a continuum/overlap (= co-expression)?
p_coexp <- ggplot(clust8@meta.data,
                  aes(x = cTEC_score, y = mTEC_I_score, color = subcluster)) +
  geom_point(size = 0.8, alpha = 0.6) +
  labs(title = "cTEC vs mTEC I module score, per cell (cluster 8 only)",
       x = "cTEC score (PSMB11/PRSS16/LY75)",
       y = "mTEC I score (KRT14/CCL19)") +
  theme_minimal()
print(p_coexp)
ggsave(paste0(DATE_PREFIX, "-cluster8-cTEC-vs-mTECI-coexpression.pdf"), p_coexp,
       width = 7, height = 5)

# Violin plots of both scores by sub-cluster -- direct, simple comparison
p_vln_cTEC <- VlnPlot_scCustom(clust8, features = "cTEC_score",
                               group.by = "subcluster", pt.size = 0.1) +
  labs(title = "cTEC marker score by sub-cluster")
p_vln_mTEC <- VlnPlot_scCustom(clust8, features = "mTEC_I_score",
                               group.by = "subcluster", pt.size = 0.1) +
  labs(title = "mTEC I marker score by sub-cluster")
p_vln_combo <- p_vln_cTEC | p_vln_mTEC
print(p_vln_combo)
ggsave(paste0(DATE_PREFIX, "-cluster8-marker-scores-violin.pdf"), p_vln_combo,
       width = 10, height = 5)

# Individual gene dotplot -- shows raw marker expression, not just composite
p_dot <- DotPlot_scCustom(clust8, features = c(cTEC_markers, mTEC_I_markers),
                          group.by = "subcluster") +
  RotatedAxis() +
  labs(title = "Individual cTEC + mTEC I markers by sub-cluster")
print(p_dot)
ggsave(paste0(DATE_PREFIX, "-cluster8-markers-dotplot.pdf"), p_dot,
       width = 6, height = 4)

# -----------------------------------------------------------------------------
# 6. Quantify: what fraction of cells co-express BOTH marker sets at a
#    meaningful level? (vs. cleanly one or the other)
# -----------------------------------------------------------------------------
# Threshold: score > 0 as "positive" for that program (module scores are
# roughly centered at 0 by construction)
clust8$cTEC_pos   <- clust8$cTEC_score > 0
clust8$mTEC_I_pos <- clust8$mTEC_I_score > 0

coexp_table <- clust8@meta.data %>%
  count(subcluster, cTEC_pos, mTEC_I_pos) %>%
  group_by(subcluster) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  ungroup()

cat("\n=== Co-expression breakdown per sub-cluster ===\n")
print(coexp_table)

write.csv(coexp_table, paste0(DATE_PREFIX, "-cluster8-coexpression-breakdown.csv"),
          row.names = FALSE)

cat("\nINTERPRETATION:\n")
cat("  High % of cTEC_pos=TRUE & mTEC_I_pos=TRUE together (double-positive)\n")
cat("  = genuine co-expression, supporting a transitional/intermediate\n")
cat("  identity rather than clean cTEC-vs-mTEC discrimination via CD205.\n")
cat("  If the cTEC-adjacent sub-cluster is mostly cTEC_pos-only (not double-\n")
cat("  positive), that argues more for a resolution/annotation artifact\n")
cat("  (misclassified cTECs) than true co-expression.\n")

# -----------------------------------------------------------------------------
# 7. Save
# -----------------------------------------------------------------------------
saveRDS(clust8, paste0(DATE_PREFIX, "-cluster8-subclustered.rds"))
