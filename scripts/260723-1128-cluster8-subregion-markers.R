# =============================================================================
# Cluster 8: gate into 4 spatial sub-regions FIRST, then re-run marker
# analysis PER sub-region -- testing whether the "mixture" signature from
# the pooled cluster-8 FindMarkers was an averaging artifact across four
# distinct populations.
# =============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(scCustomize)

DATE_PREFIX <- format(Sys.time(), "%y%m%d-%H%M")

tec <- readRDS("data/rds-objects/260501-after-annotation-1.rds")
Idents(tec) <- "seurat_clusters"

# -----------------------------------------------------------------------------
# STEP 1: Coordinate-based gating of cluster 8 into 4 sub-regions
# -----------------------------------------------------------------------------
umap_coords <- Embeddings(tec, "umap") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell") %>%
  rename(umap_1 = 2, umap_2 = 3)   # ensure consistent column names regardless
# of Seurat's default UMAP_1/UMAP_2 naming
umap_coords$seurat_clusters <- as.character(tec$seurat_clusters[umap_coords$cell])

c8_coords <- umap_coords %>% filter(seurat_clusters == "8")

# ADJUST these bounding boxes based on visual inspection of your c8-only
# aggregate DimPlot (approximate values from what you shared):
c8_coords <- c8_coords %>%
  mutate(sub_region = case_when(
    umap_1 >= -8  & umap_1 <=  0 & umap_2 >=  1  & umap_2 <=  5  ~ "near_cTEC",
    umap_1 >=  5  & umap_1 <= 11 & umap_2 >=  3  & umap_2 <=  9  ~ "near_Nurse",
    umap_1 >= -5  & umap_1 <=  2 & umap_2 >= -7.5  & umap_2 <= -3  ~ "near_mTEC_I",
    umap_1 >= -3  & umap_1 <=  5 & umap_2 >= -13 & umap_2 <= -8  ~ "near_mTEC_II",
    TRUE ~ "unassigned"
  ))

cat("=== Cluster 8 sub-region cell counts ===\n")
print(table(c8_coords$sub_region))

tec$c8_subregion <- NA_character_
tec$c8_subregion[match(c8_coords$cell, colnames(tec))] <- c8_coords$sub_region

# Visual check -- confirm the gating captured what you intended
p_gating <- DimPlot_scCustom(subset(tec, subset = seurat_clusters == "8"),
                             group.by = "c8_subregion", label = TRUE) +
  labs(title = "Cluster 8 gated into 4 sub-regions")
print(p_gating)
ggsave(paste0(DATE_PREFIX, "-cluster8-gating-check.pdf"), p_gating,
       width = 8, height = 6)

# -----------------------------------------------------------------------------
# STEP 2: Define neighbor cluster mapping for the one-vs-specific-neighbor test
# -----------------------------------------------------------------------------
# ADJUST -- numeric pre-annotation cluster IDs for each neighbor cell type.
# Check your new.cluster.ids mapping vector to confirm these.
neighbor_map <- list(
  near_cTEC    = c("0", "1"),
  near_mTEC_I  = c("2"),
  near_mTEC_II = c("7", "18"),
  near_Nurse   = c("3", "14")
)

# -----------------------------------------------------------------------------
# STEP 3: Per sub-region -- create a temporary identity column, then run
# both marker analyses
# -----------------------------------------------------------------------------
# Give every cell a per-subregion identity: cluster-8 cells get their
# sub-region label, everyone else keeps their normal cluster ID. This lets
# FindMarkers treat each sub-region as its own group correctly.
tec$ident_for_c8_analysis <- ifelse(
  !is.na(tec$c8_subregion),
  tec$c8_subregion,
  as.character(tec$seurat_clusters)
)
Idents(tec) <- "ident_for_c8_analysis"

results_all       <- list()
results_neighbor  <- list()

for (region in names(neighbor_map)) {
  cat("\n", strrep("=", 60), "\n")
  cat("SUB-REGION:", region, "\n")
  cat(strrep("=", 60), "\n")
  
  n_cells <- sum(tec$ident_for_c8_analysis == region, na.rm = TRUE)
  cat("Cell count:", n_cells, "\n")
  if (n_cells < 20) {
    cat("WARNING: very few cells -- results will be noisy/low-powered.\n")
  }
  
  # --- One-vs-all (this sub-region vs. everything else in the dataset) ---
  m_all <- tryCatch({
    FindMarkers(
      tec, ident.1 = region, only.pos = TRUE,
      min.pct = 0.25, logfc.threshold = 0.25, verbose = TRUE
    )
  }, error = function(e) { cat("FAILED (likely too few cells):", conditionMessage(e), "\n"); NULL })
  
  if (!is.null(m_all)) {
    m_all$gene <- rownames(m_all)
    m_all <- m_all %>% arrange(desc(avg_log2FC))
    results_all[[region]] <- m_all
    cat("\nTop 15, vs-ALL:\n")
    print(head(m_all[, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")], 15))
    write.csv(m_all, paste0(DATE_PREFIX, "-c8-", region, "-vs-ALL.csv"),
              row.names = FALSE)
  }
  
  # --- One-vs-specific-neighbor (this sub-region vs. ONLY its own neighbor) ---
  m_neighbor <- tryCatch({
    FindMarkers(
      tec, ident.1 = region, ident.2 = neighbor_map[[region]],
      only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = TRUE
    )
  }, error = function(e) { cat("FAILED (likely too few cells):", conditionMessage(e), "\n"); NULL })
  
  if (!is.null(m_neighbor)) {
    m_neighbor$gene <- rownames(m_neighbor)
    m_neighbor <- m_neighbor %>% arrange(desc(avg_log2FC))
    results_neighbor[[region]] <- m_neighbor
    cat("\nTop 15, vs-SPECIFIC-NEIGHBOR (", paste(neighbor_map[[region]], collapse=","), "):\n")
    print(head(m_neighbor[, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")], 15))
    write.csv(m_neighbor, paste0(DATE_PREFIX, "-c8-", region, "-vs-NEIGHBOR.csv"),
              row.names = FALSE)
  }
}

# -----------------------------------------------------------------------------
# STEP 4: Cross-sub-region comparison -- did the "mixture" genes concentrate
# in one specific sub-region?
# -----------------------------------------------------------------------------
mixture_genes <- c("POU2AF1", "INSM1", "STAR", "SV2B", "CALML3", "FOXQ1",
                   "FOXS1", "TRAF1", "RELT", "LTB", "IL21R", "FXYD2", "WNT4")

cat("\n", strrep("=", 60), "\n")
cat("Where did the original 'mixture' genes end up? (top-30 rank per sub-region)\n")
cat(strrep("=", 60), "\n")
for (region in names(results_all)) {
  ranks <- match(mixture_genes, results_all[[region]]$gene)
  names(ranks) <- mixture_genes
  cat("\n", region, ":\n")
  print(ranks[!is.na(ranks)])
}