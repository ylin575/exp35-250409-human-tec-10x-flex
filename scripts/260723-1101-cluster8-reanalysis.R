# =============================================================================
# Cluster 8 re-analysis -- FOUR PHASES + preliminary marker-identification
#
# Motivation: cluster 8 has FOUR spatially-distinct sub-regions on the joint
# UMAP, each near a different neighbor (cTEC, mTEC I, mTEC II, Nurse), not
# two as originally assumed. Also, we've never actually asked what genes
# make cluster 8 stand out on its own -- that comes first.
# =============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(scCustomize)
library(harmony)

DATE_PREFIX <- format(Sys.time(), "%y%m%d-%H%M")

tec <- readRDS("data/rds-objects/260501-after-annotation-1.rds")

# -----------------------------------------------------------------------------
# PHASE 0 -- What defines cluster 8?  Two complementary marker analyses.
# -----------------------------------------------------------------------------
# Assumes seurat_clusters holds the 27 pre-annotation Leiden IDs.
# Confirm before running:
stopifnot("seurat_clusters" %in% colnames(tec@meta.data))
Idents(tec) <- "seurat_clusters"

# --- 0a. One-vs-all: cluster 8 vs. everything else combined ---
# Answers: "what genes are elevated in cluster 8 relative to the whole
# dataset averaged together?"  Cluster 8's distinctive positive markers.
markers_c8_vs_all <- FindMarkers(
  tec,
  ident.1         = "8",
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25,
  verbose         = TRUE
)
markers_c8_vs_all$gene <- rownames(markers_c8_vs_all)
markers_c8_vs_all <- markers_c8_vs_all %>%
  arrange(desc(avg_log2FC))

write.csv(markers_c8_vs_all,
          paste0(DATE_PREFIX, "-cluster8-markers-vs-ALL.csv"),
          row.names = FALSE)

cat("=== Top 30 genes: cluster 8 vs. all other clusters ===\n")
print(head(markers_c8_vs_all[, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")], 30))

# --- 0b. Pairwise: cluster 8 vs. its four spatial neighbors ---
# Answers: "what specifically separates cluster 8 from cTEC, mTEC I, mTEC II,
# and Nurse -- the populations sitting adjacent to it in UMAP space?"
# This is closer to what actually makes cluster 8 a distinct cluster given
# its spatial context. Adjust ident.2 based on which numeric clusters
# correspond to those four cell types in YOUR annotation mapping.

# ADJUST these -- these are the NUMERIC pre-annotation cluster IDs for
# cTEC, mTEC I, mTEC II, and Nurse. Look at new.cluster.ids in your
# original clustering script to confirm.
NEIGHBOR_CLUSTERS <- c("0", "1", "22",   # <-- adjust: cTEC-annotated clusters
                       "6",              # <-- adjust: mTEC I proper (excluding 8)
                       "7", "18",        # <-- adjust: mTEC II
                       "3", "14")        # <-- adjust: Nurse

markers_c8_vs_neighbors <- FindMarkers(
  tec,
  ident.1         = "8",
  ident.2         = NEIGHBOR_CLUSTERS,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25,
  verbose         = TRUE
)
markers_c8_vs_neighbors$gene <- rownames(markers_c8_vs_neighbors)
markers_c8_vs_neighbors <- markers_c8_vs_neighbors %>%
  arrange(desc(avg_log2FC))

write.csv(markers_c8_vs_neighbors,
          paste0(DATE_PREFIX, "-cluster8-markers-vs-NEIGHBORS.csv"),
          row.names = FALSE)

cat("\n=== Top 30 genes: cluster 8 vs. spatial neighbors only ===\n")
print(head(markers_c8_vs_neighbors[, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")], 30))

# --- 0c. Overlap analysis: shared vs. unique markers between the two views ---
top_all       <- head(markers_c8_vs_all$gene, 30)
top_neighbors <- head(markers_c8_vs_neighbors$gene, 30)

cat("\n=== Markers common to both views (cluster 8's core signature) ===\n")
print(intersect(top_all, top_neighbors))
cat("\n=== Markers unique to vs-ALL (broad character, less discriminating) ===\n")
print(setdiff(top_all, top_neighbors))
cat("\n=== Markers unique to vs-NEIGHBORS (sibling-discriminating) ===\n")
print(setdiff(top_neighbors, top_all))

# --- 0d. Build a "cluster 8 identity" module score from these markers ---
# Take the top N genes from the pairwise analysis (best discriminators
# from spatial neighbors, which is the analysis-relevant view here)
c8_module_genes <- head(markers_c8_vs_neighbors$gene, 15)

tec <- AddModuleScore(
  tec,
  features = list(c8_module_genes),
  name     = "C8_identity",
  verbose  = TRUE
)
tec$C8_identity_score <- tec$C8_identity1
tec$C8_identity1 <- NULL

# Where does this "cluster 8 identity" signal actually concentrate?
# If it truly identifies cluster 8, it should peak in cluster 8 specifically.
Idents(tec) <- "seurat_clusters"
p_c8mod <- VlnPlot_scCustom(tec, features = "C8_identity_score",
                            group.by = "seurat_clusters", pt.size = 0) +
  labs(title = "Cluster 8 identity module score across all clusters")
print(p_c8mod)
ggsave(paste0(DATE_PREFIX, "-cluster8-identity-module-score-vlnplot.pdf"),
       p_c8mod, width = 14, height = 5)


# -----------------------------------------------------------------------------
# PHASE 1 -- Coordinate-based bounding box gating of cluster 8's four regions
# -----------------------------------------------------------------------------
# Approximate coordinates from your UMAPs -- ADJUST after visual inspection.

umap_coords <- Embeddings(tec, "umap") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell") %>%
  mutate(seurat_clusters = as.character(tec$seurat_clusters[cell]))

c8_coords <- umap_coords %>% filter(seurat_clusters == "8")

# ADJUST these bounding boxes based on visual inspection of your c8-only
# DimPlot from the last upload. Approximate values from that plot:
c8_coords <- c8_coords %>%
  mutate(sub_region = case_when(
    umap_1 >= -8  & umap_1 <=  0 & umap_2 >=  0 & umap_2 <=  6  ~ "near_cTEC",       # upper-left
    umap_1 >=  5  & umap_1 <= 10 & umap_2 >=  5 & umap_2 <=  9  ~ "near_Nurse",       # upper-right
    umap_1 >= -3  & umap_1 <=  0 & umap_2 >= -6 & umap_2 <= -3  ~ "near_mTEC_I",     # middle core
    umap_1 >= -3  & umap_1 <=  3 & umap_2 >= -12 & umap_2 <= -8 ~ "near_mTEC_II",    # lower blob
    TRUE ~ "unassigned"
  ))

cat("\n=== Cluster 8 sub-region assignments ===\n")
print(table(c8_coords$sub_region))

# Attach back to tec meta.data
tec$c8_subregion <- NA_character_
tec$c8_subregion[match(c8_coords$cell, colnames(tec))] <- c8_coords$sub_region

# Visualize the gating result to confirm the boxes captured what you intended
p_gating <- DimPlot_scCustom(subset(tec, subset = seurat_clusters == "8"),
                             group.by = "c8_subregion", label = TRUE) +
  labs(title = "Cluster 8 sub-regions after coordinate-based gating")
print(p_gating)
ggsave(paste0(DATE_PREFIX, "-cluster8-subregion-gating.pdf"),
       p_gating, width = 8, height = 6)


# -----------------------------------------------------------------------------
# PHASE 2 -- Four module scores based on cell_type annotation
# -----------------------------------------------------------------------------
# Marker panels for each neighbor identity. Adjust based on YOUR manuscript's
# established markers for each cell type.

cTEC_markers    <- c("PSMB11", "PRSS16", "LY75")
mTEC_I_markers  <- c("KRT14", "CCL19", "CCL21")   # KRT5 added optionally below
mTEC_II_markers <- c("AIRE", "FEZF2", "SPIB")
Nurse_markers   <- c("CD52", "PTPRC", "CORO1A")   # ADJUST -- use YOUR nurse markers

# Filter to genes actually present in the object
cTEC_markers    <- intersect(cTEC_markers,    rownames(tec))
mTEC_I_markers  <- intersect(mTEC_I_markers,  rownames(tec))
mTEC_II_markers <- intersect(mTEC_II_markers, rownames(tec))
Nurse_markers   <- intersect(Nurse_markers,   rownames(tec))

for (mod in list(
  list(name = "cTEC_score",    genes = cTEC_markers),
  list(name = "mTEC_I_score",  genes = mTEC_I_markers),
  list(name = "mTEC_II_score", genes = mTEC_II_markers),
  list(name = "Nurse_score",   genes = Nurse_markers)
)) {
  tec <- AddModuleScore(tec, features = list(mod$genes),
                        name = mod$name, verbose = TRUE)
  tec[[mod$name]] <- tec[[paste0(mod$name, "1")]]
  tec[[paste0(mod$name, "1")]] <- NULL
}


# -----------------------------------------------------------------------------
# PHASE 3 -- Per-sub-region co-expression analysis
# -----------------------------------------------------------------------------
# For each sub-region: what's its own identity score, its spatial neighbor's
# score, and the co-expression rate (soft + strict)?

subregion_summary <- tec@meta.data %>%
  filter(!is.na(c8_subregion)) %>%
  group_by(c8_subregion) %>%
  summarise(
    n = n(),
    mean_cTEC    = round(mean(cTEC_score),    3),
    mean_mTEC_I  = round(mean(mTEC_I_score),  3),
    mean_mTEC_II = round(mean(mTEC_II_score), 3),
    mean_Nurse   = round(mean(Nurse_score),   3),
    .groups = "drop"
  )

cat("\n=== Mean module scores per cluster-8 sub-region ===\n")
print(subregion_summary)

write.csv(subregion_summary,
          paste0(DATE_PREFIX, "-cluster8-subregion-module-scores.csv"),
          row.names = FALSE)

# The reviewer's specific concern: near_cTEC sub-region co-expression
# using the SAME strict raw-UMI criterion as before for consistency
counts_mat <- GetAssayData(tec, assay = "RNA", layer = "counts")
MIN_UMI <- 2

for (subregion in unique(na.omit(tec$c8_subregion))) {
  cells <- colnames(tec)[!is.na(tec$c8_subregion) & tec$c8_subregion == subregion]
  cTEC_sum   <- Matrix::colSums(counts_mat[cTEC_markers,   cells, drop = FALSE])
  mTEC_I_sum <- Matrix::colSums(counts_mat[mTEC_I_markers, cells, drop = FALSE])
  double_pos <- sum(cTEC_sum >= MIN_UMI & mTEC_I_sum >= MIN_UMI)
  total <- length(cells)
  cat(sprintf("Sub-region %s: %d/%d (%.1f%%) cTEC+/mTEC-I+ double-positive (raw UMI >= %d)\n",
              subregion, double_pos, total, 100 * double_pos / total, MIN_UMI))
}


# -----------------------------------------------------------------------------
# PHASE 4 (INDEPENDENT) -- Harmony integration of aggregate sample
# -----------------------------------------------------------------------------
# Tag the current cluster-8 cells so we can track where they land after
# integration + re-clustering. This is the whole point of the phase --
# does cluster 8 hold together as its own cluster post-integration?

tec$original_c8_flag <- ifelse(tec$seurat_clusters == "8",
                               "cluster_8", "other")

DefaultAssay(tec) <- "RNA"

tec_harmony <- tec %>%
  NormalizeData(verbose = TRUE) %>%
  FindVariableFeatures(nfeatures = 2000, verbose = TRUE) %>%
  ScaleData(verbose = TRUE) %>%
  RunPCA(npcs = 30, verbose = TRUE) %>%
  RunHarmony(
    group.by.vars  = "sample.id",
    reduction      = "pca",
    dims.use       = 1:30,
    reduction.save = "harmony",
    verbose        = TRUE
  ) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20, verbose = TRUE) %>%
  FindClusters(resolution = 0.8, verbose = TRUE) %>%
  RunUMAP(reduction = "harmony", dims = 1:20,
          reduction.name = "umap_harmony", verbose = TRUE)

# Where do the original cluster-8 cells land post-integration?
cat("\n=== Post-Harmony cluster distribution of original cluster-8 cells ===\n")
print(table(tec_harmony$seurat_clusters[tec_harmony$original_c8_flag == "cluster_8"]))

# Visualize
p_harmony <- DimPlot_scCustom(tec_harmony, reduction = "umap_harmony",
                              group.by = "original_c8_flag") +
  labs(title = "Post-Harmony UMAP -- original cluster 8 cells highlighted")
print(p_harmony)
ggsave(paste0(DATE_PREFIX, "-cluster8-post-harmony-tracking.pdf"),
       p_harmony, width = 8, height = 6)

# Are the double-positive cells specifically still co-clustered, or scattered?
p_harmony_subregions <- DimPlot_scCustom(
  tec_harmony, reduction = "umap_harmony",
  group.by = "c8_subregion", na.value = "grey90"
) +
  labs(title = "Post-Harmony UMAP -- cluster 8 sub-regions")
print(p_harmony_subregions)
ggsave(paste0(DATE_PREFIX, "-cluster8-subregions-post-harmony.pdf"),
       p_harmony_subregions, width = 8, height = 6)

saveRDS(tec_harmony,
        paste0(DATE_PREFIX, "-tec-harmony-cluster8-tracked.rds"))