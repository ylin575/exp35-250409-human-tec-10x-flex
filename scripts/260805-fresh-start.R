# =============================================================================
# Combined unintegrated + Harmony-integrated pipeline, built fresh from raw
# counts, with doublet detection/removal included from the start.
# =============================================================================

library(Seurat)
library(harmony)
library(scDblFinder)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(patchwork)

DATE_PREFIX <- format(Sys.time(), "%y%m%d-%H%M")
set.seed(42)





# -----------------------------------------------------------------------------
# Load counts (adjust paths to match your directory structure)
# -----------------------------------------------------------------------------
ht67_cd205neg.data <- Read10X_h5("rawdata/h5-files/ht67-cd205neg/sample_filtered_feature_bc_matrix.h5")
ht67_cd205pos.data <- Read10X_h5("rawdata/h5-files/ht67-cd205pos/sample_filtered_feature_bc_matrix.h5")
ht70.data          <- Read10X_h5("rawdata/h5-files/ht70/sample_filtered_feature_bc_matrix.h5")
ht71.data          <- Read10X_h5("rawdata/h5-files/ht71/sample_filtered_feature_bc_matrix.h5")





# -----------------------------------------------------------------------------
# Create Seurat objects
# -----------------------------------------------------------------------------
ht67_cd205neg <- CreateSeuratObject(counts = ht67_cd205neg.data, project = "htec-10xflex",
                                    min.cells = 3, min.features = 200)
ht67_cd205pos <- CreateSeuratObject(counts = ht67_cd205pos.data, project = "htec-10xflex",
                                    min.cells = 3, min.features = 200)
ht70 <- CreateSeuratObject(counts = ht70.data, project = "htec-10xflex",
                           min.cells = 3, min.features = 200)
ht71 <- CreateSeuratObject(counts = ht71.data, project = "htec-10xflex",
                           min.cells = 3, min.features = 200)

# sample.id and donor metadata
ht67_cd205neg$sample.id <- "ht67-cd205neg"; ht67_cd205neg$donor <- "HT67"
ht67_cd205pos$sample.id <- "ht67-cd205pos"; ht67_cd205pos$donor <- "HT67"
ht70$sample.id <- "ht70"; ht70$donor <- "HT70"
ht71$sample.id <- "ht71"; ht71$donor <- "HT71"

# percent.mt
ht67_cd205neg[["percent.mt"]] <- PercentageFeatureSet(ht67_cd205neg, pattern = "^MT-")
ht67_cd205pos[["percent.mt"]] <- PercentageFeatureSet(ht67_cd205pos, pattern = "^MT-")
ht70[["percent.mt"]] <- PercentageFeatureSet(ht70, pattern = "^MT-")
ht71[["percent.mt"]] <- PercentageFeatureSet(ht71, pattern = "^MT-")





# -----------------------------------------------------------------------------
# Per-sample QC filtering (mean + 3*SD, matching original script)
# -----------------------------------------------------------------------------
qc_filter <- function(obj) {
  upper_mt <- mean(obj$percent.mt) + 3 * sd(obj$percent.mt)
  upper_feature <- mean(obj$nFeature_RNA) + 3 * sd(obj$nFeature_RNA)
  obj <- subset(obj, subset = percent.mt <= upper_mt)
  obj <- subset(obj, subset = nFeature_RNA >= 500 & nFeature_RNA <= upper_feature)
  obj
}

ht67_cd205neg_sub <- qc_filter(ht67_cd205neg)
ht67_cd205pos_sub <- qc_filter(ht67_cd205pos)
ht70_sub           <- qc_filter(ht70)
ht71_sub           <- qc_filter(ht71)





# -----------------------------------------------------------------------------
# Merge, join layers
# -----------------------------------------------------------------------------
samples <- merge(x = ht67_cd205neg_sub,
                 y = list(ht67_cd205pos_sub, ht70_sub, ht71_sub))
samples <- JoinLayers(samples)

table(samples$sample.id, samples$donor)   # sanity check





# -----------------------------------------------------------------------------
# Standard preprocessing (shared by both branches)
# -----------------------------------------------------------------------------
samples <- NormalizeData(samples)
samples <- FindVariableFeatures(samples)
samples <- ScaleData(samples)
samples <- RunPCA(samples)

ElbowPlot(samples, ndims = 60)





# -----------------------------------------------------------------------------
# BRANCH A: Unintegrated
# -----------------------------------------------------------------------------
samples <- FindNeighbors(samples, reduction = "pca", dims = 1:18,
                         graph.name = c("nn_unintegrated", "snn_unintegrated"))
samples <- FindClusters(samples, graph.name = "snn_unintegrated",
                        resolution = 0.5, cluster.name = "clusters_unintegrated")
samples <- RunUMAP(samples, reduction = "pca", dims = 1:18,
                   reduction.name = "umap_unintegrated")





# -----------------------------------------------------------------------------
# BRANCH B: Harmony-integrated
# -----------------------------------------------------------------------------
samples <- RunHarmony(samples, group.by.vars = "donor", reduction = "pca",
                      dims.use = 1:18, reduction.save = "harmony")
samples <- FindNeighbors(samples, reduction = "harmony", dims = 1:18,
                         graph.name = c("nn_harmony", "snn_harmony"))
samples <- FindClusters(samples, graph.name = "snn_harmony",
                        resolution = 0.5, cluster.name = "clusters_harmony")
samples <- RunUMAP(samples, reduction = "harmony", dims = 1:18,
                   reduction.name = "umap_harmony")





# -----------------------------------------------------------------------------
# Sanity checks
# -----------------------------------------------------------------------------
cat("Reductions present:", paste(Reductions(samples), collapse = ", "), "\n")
cat("Unintegrated clusters found:", length(unique(samples$clusters_unintegrated)), "\n")
cat("Harmony clusters found:", length(unique(samples$clusters_harmony)), "\n")
stopifnot(all(c("umap_unintegrated", "umap_harmony", "harmony") %in% Reductions(samples)))
stopifnot(!any(is.na(samples$clusters_unintegrated)))
stopifnot(!any(is.na(samples$clusters_harmony)))

table(samples$clusters_unintegrated, samples$clusters_harmony)

#saveRDS(samples, paste0(DATE_PREFIX, "-samples-both-pre-and-post-integration.rds"))
#samples <- readRDS(file = "data/rds-objects/260805-1608-samples-both-pre-and-post-integration.rds")





# -----------------------------------------------------------------------------
# Find markers
# -----------------------------------------------------------------------------
Idents(samples) <- "clusters_unintegrated"
stopifnot(identical(as.character(Idents(samples)), as.character(samples$clusters_unintegrated)))
markers_unint <- FindAllMarkers(samples, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

Idents(samples) <- "clusters_harmony"
stopifnot(identical(as.character(Idents(samples)), as.character(samples$clusters_harmony)))
markers_harmony <- FindAllMarkers(samples, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


# Full marker tables
write.csv(markers_unint, paste0(DATE_PREFIX, "-unintegrated-markers-FULL.csv"), row.names = FALSE)
write.csv(markers_harmony, paste0(DATE_PREFIX, "-harmony-markers-FULL.csv"), row.names = FALSE)

# Top 15 per cluster, for a quick scan
top15_unint <- markers_unint %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 15) %>%
  ungroup()

top15_harmony <- markers_harmony %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 15) %>%
  ungroup()

write.csv(top15_unint, paste0(DATE_PREFIX, "-unintegrated-markers-TOP15.csv"), row.names = FALSE)
write.csv(top15_harmony, paste0(DATE_PREFIX, "-harmony-markers-TOP15.csv"), row.names = FALSE)

#saveRDS(samples, paste0(DATE_PREFIX, "-samples-both-pre-and-post-integration-post-find-markers.rds"))
samples <- readRDS(file = "data/rds-objects/260805-2223-samples-both-pre-and-post-integration-post-find-markers.rds")





# -----------------------------------------------------------------------------
# Annotating cell types
# -----------------------------------------------------------------------------
set.seed(2)

# annotations for unintegrated data
Idents(samples) <- "clusters_unintegrated"

new.cluster.ids <- c(
  "cTEC",                         
  "mTEC lo",                        
  "cTEC",                         
  "Nurse",                         
  "BMC TEC",                      
  "BMC TEC",                   
  "mTEC hi",                 
  "Bipolar TEC",                  
  "mTEC lo",               
  "BMC TEC",                      
  "Mimetic",                      
  "Mimetic",                      
  "Mimetic",                      
  "Nurse",                         
  "Mimetic",                      
  "Nurse",                        
  "Pericyte",                    
  "Mimetic",                     
  "Mimetic",                      
  "Mimetic",                       
  "Endothelial",                  
  "Fibroblast",                   
  "Erythroid",                   
  "Mimetic",                     
  "Erythroid",                   
  "Mimetic"                   
)

# Sanity check: vector length must exactly match number of clusters, in order
stopifnot(length(new.cluster.ids) == length(levels(samples)))
cat("Cluster levels (order matters!):\n")
print(levels(samples))
cat("\nProposed annotation, aligned:\n")
print(data.frame(cluster = levels(samples), annotation = new.cluster.ids))

names(new.cluster.ids) <- levels(samples)
samples <- RenameIdents(samples, new.cluster.ids)
samples$cell_type_unintegrated <- Idents(samples)

# Restore numeric cluster identity as the active Idents (so downstream code
# that expects cluster numbers, e.g. FindAllMarkers re-runs, still works)
Idents(samples) <- "clusters_unintegrated"

# Final sanity check
stopifnot(!any(is.na(samples$cell_type_unintegrated)))
cat("\n=== Final annotation counts ===\n")
print(table(samples$cell_type_unintegrated))

# check annotation on umap
DimPlot(samples, reduction = "umap_unintegrated",
        group.by = "cell_type_unintegrated", label = TRUE, repel = TRUE) +
  labs(title = "Unintegrated: annotated cell types")

# annotations for integrated data
Idents(samples) <- "clusters_harmony"

new.cluster.ids.harmony <- c(
  "cTEC",                     
  "cTEC",                         
  "mTEC lo",                      
  "BMC TEC",                     
  "Nurse",                        
  "mTEC lo",                      
  "Bipolar TEC",                  
  "mTEC hi",                     
  "mTEC lo",                     
  "Mimetic",                     
  "Mimetic",                     
  "Nurse",                       
  "Mimetic",                    
  "Mimetic",                     
  "Mimetic",                      
  "Pericyte",                  
  "Nurse",                      
  "Erythroid",                   
  "Mimetic",                      
  "mTEC hi", 
  "Endothelial",                 
  "Mimetic"                      
)

# Sanity check: vector length must exactly match number of clusters, in order
stopifnot(length(new.cluster.ids.harmony) == length(levels(samples)))
cat("Cluster levels (order matters!):\n")
print(levels(samples))
cat("\nProposed annotation, aligned:\n")
print(data.frame(cluster = levels(samples), annotation = new.cluster.ids.harmony))

names(new.cluster.ids.harmony) <- levels(samples)
samples <- RenameIdents(samples, new.cluster.ids.harmony)
samples$cell_type_harmony <- Idents(samples)

# Restore numeric cluster identity as active Idents
Idents(samples) <- "clusters_harmony"

# Final sanity check
stopifnot(!any(is.na(samples$cell_type_harmony)))
cat("\n=== Final annotation counts ===\n")
print(table(samples$cell_type_harmony))

# check annotation on umap
DimPlot(samples, reduction = "umap_harmony",
        group.by = "cell_type_harmony", label = TRUE, repel = TRUE) +
  labs(title = "Harmony-integrated: annotated cell types")

#saveRDS(samples, paste0(DATE_PREFIX, "-samples-both-pre-and-post-integration-post-annotating.rds"))
#samples <- readRDS(file = "data/rds-objects/260806-1121-samples-both-pre-and-post-integration-post-annotating.rds")





# -----------------------------------------------------------------------------
# Ragazzini PolyKRT module score -- unintegrated data
# -----------------------------------------------------------------------------
polyKRT_module_genes <- c(
  "CEBPD", "CLU", "FN1", "IFITM3", "TIMP1", "VCAM1", "TAGLN",
  "BCAM", "LIFR", "CH25H", "CCN2", "CCL19",
  "KRT13", "KRT14", "KRT15", "KRT17", "KRT19"
)

missing_genes <- setdiff(polyKRT_module_genes, rownames(samples))
if (length(missing_genes) > 0) {
  warning("Not in object, dropping: ", paste(missing_genes, collapse = ", "))
  polyKRT_module_genes <- setdiff(polyKRT_module_genes, missing_genes)
}
cat("Module genes used (n=", length(polyKRT_module_genes), "):\n", sep = "")
cat(paste(polyKRT_module_genes, collapse = ", "), "\n")

samples <- AddModuleScore(
  samples, features = list(polyKRT_module_genes),
  name = "PolyKRT_module", verbose = TRUE
)
samples$PolyKRT_module_score <- samples$PolyKRT_module1
samples$PolyKRT_module1 <- NULL

# Sanity check
stopifnot(!any(is.na(samples$PolyKRT_module_score)))
cat("Score range:", round(range(samples$PolyKRT_module_score), 3), "\n")





# -----------------------------------------------------------------------------
# Cell type order (adjust if you prefer a different arrangement)
# -----------------------------------------------------------------------------
cluster_order_unint <- c("Nurse", "cTEC", "BMC TEC", "mTEC lo", "mTEC hi",
                         "Mimetic", "Endothelial","Bipolar TEC","Pericyte",
                         "Fibroblast","Erythroid")

stopifnot(setequal(cluster_order_unint, levels(samples$cell_type_unintegrated)))
samples$cell_type_unintegrated <- factor(samples$cell_type_unintegrated,
                                         levels = cluster_order_unint)




# -----------------------------------------------------------------------------
# VlnPlot with a subset of cell types excluded
# -----------------------------------------------------------------------------
exclude_types <- c("Bipolar TEC", "Pericyte", "Fibroblast", "Erythroid")

samples_vln_subset <- subset(
  samples,
  subset = !(cell_type_unintegrated %in% exclude_types)
)

# Drop unused factor levels so they don't show as empty categories on the axis
samples_vln_subset$cell_type_unintegrated <- droplevels(samples_vln_subset$cell_type_unintegrated)

# Sanity check
cat("Cell types remaining:\n")
print(table(samples_vln_subset$cell_type_unintegrated))
stopifnot(!any(levels(samples_vln_subset$cell_type_unintegrated) %in% exclude_types))

p_vln_subset <- VlnPlot(samples_vln_subset, features = "PolyKRT_module_score",
                        group.by = "cell_type_unintegrated", pt.size = 0.1) +
  labs(title = "Ragazzini PolyKRT module score by cell type (unintegrated, subset)",
       y = "Module score")
print(p_vln_subset)
ggsave(paste0(DATE_PREFIX, "-PolyKRT-module-score-violin-unintegrated-subset.pdf"),
       p_vln_subset, width = 8, height = 5)





# -----------------------------------------------------------------------------
# 2. FeaturePlot
# -----------------------------------------------------------------------------
p_umap <- FeaturePlot(samples, features = "PolyKRT_module_score",
                      reduction = "umap_unintegrated", order = TRUE) +
  scale_color_viridis_c(option = "viridis") +
  labs(title = "PolyKRT module score, UMAP (unintegrated)")
print(p_umap)
ggsave(paste0(DATE_PREFIX, "-PolyKRT-module-score-umap-unintegrated.pdf"),
       p_umap, width = 6, height = 5)





# -----------------------------------------------------------------------------
# Heatmap, subset (Bipolar TEC, Pericyte, Fibroblast, Erythroid excluded)
# -----------------------------------------------------------------------------
cluster_order_subset <- c("Nurse", "cTEC", "BMC TEC", "mTEC lo", "mTEC hi",
                          "Mimetic", "Endothelial")

stopifnot(setequal(cluster_order_subset,
                   levels(samples_vln_subset$cell_type_unintegrated)))

avg_mat_subset <- AverageExpression(
  samples_vln_subset, features = polyKRT_module_genes,
  assay = "RNA", group.by = "cell_type_unintegrated"
)$RNA
avg_mat_subset <- avg_mat_subset[polyKRT_module_genes, cluster_order_subset, drop = FALSE]

log_mat_subset <- log1p(avg_mat_subset)
z_mat_subset   <- t(scale(t(log_mat_subset)))
z_mat_subset[is.nan(z_mat_subset)] <- 0

score_by_cluster_subset <- samples_vln_subset@meta.data %>%
  group_by(cell_type_unintegrated) %>%
  summarise(mean_score = mean(PolyKRT_module_score), .groups = "drop") %>%
  tibble::deframe()
score_by_cluster_subset <- score_by_cluster_subset[cluster_order_subset]

score_z_subset <- as.numeric(scale(score_by_cluster_subset))
score_z_subset[is.nan(score_z_subset)] <- 0
names(score_z_subset) <- cluster_order_subset

z_mat_full_subset <- rbind(z_mat_subset, PolyKRT_module_score = score_z_subset)
gene_order_full <- c(polyKRT_module_genes, "PolyKRT_module_score")

df_heatmap_subset <- as.data.frame(z_mat_full_subset) %>%
  tibble::rownames_to_column("gene") %>%
  tidyr::pivot_longer(-gene, names_to = "cluster", values_to = "z") %>%
  mutate(
    gene    = factor(gene, levels = rev(gene_order_full)),
    cluster = factor(cluster, levels = cluster_order_subset)
  )

p_heatmap_subset <- ggplot(df_heatmap_subset, aes(x = cluster, y = gene, fill = z)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient2(
    low = "#F7E27A", mid = "grey95", high = "#1B3B6F",
    midpoint = 0, name = "Relative\nexpr. (z)"
  ) +
  labs(title = "PolyKRT module genes + composite score, by cell type (unintegrated, subset)") +
  theme_minimal(base_size = 10) +
  theme(
    axis.title  = element_blank(),
    plot.title  = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid  = element_blank()
  )
print(p_heatmap_subset)
ggsave(paste0(DATE_PREFIX, "-PolyKRT-module-heatmap-unintegrated-subset.pdf"),
       p_heatmap_subset, width = 9, height = 6)





# -----------------------------------------------------------------------------
# Ragazzini PolyKRT Score: summary by cell type
# -----------------------------------------------------------------------------
score_summary <- samples@meta.data %>%
  group_by(cell_type_unintegrated) %>%
  summarise(
    n            = n(),
    mean_score   = round(mean(PolyKRT_module_score), 3),
    median_score = round(median(PolyKRT_module_score), 3),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_score))

cat("=== Mean PolyKRT module score per cell type (unintegrated) ===\n")
print(score_summary)

write.csv(score_summary,
          paste0(DATE_PREFIX, "-PolyKRT-module-score-by-celltype-unintegrated.csv"),
          row.names = FALSE)





# -----------------------------------------------------------------------------
# Ragazzini PolyKRT Score: Pairwise Wilcoxon tests -- top-scoring type vs every other type
# -----------------------------------------------------------------------------
top_type <- score_summary$cell_type_unintegrated[1]
other_types <- setdiff(unique(samples$cell_type_unintegrated), top_type)

pairwise_tests <- lapply(other_types, function(ct) {
  top_scores   <- samples$PolyKRT_module_score[samples$cell_type_unintegrated == top_type]
  other_scores <- samples$PolyKRT_module_score[samples$cell_type_unintegrated == ct]
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

cat("\n=== Highest-scoring cell type (", as.character(top_type),
    ") vs. every other type (Wilcoxon, BH-adjusted) ===\n", sep = "")
print(pairwise_tests)

write.csv(pairwise_tests,
          paste0(DATE_PREFIX, "-PolyKRT-module-pairwise-tests-unintegrated.csv"),
          row.names = FALSE)





# -----------------------------------------------------------------------------
# Ragazzini PolyKRT: Concentration test
# -----------------------------------------------------------------------------
N_TOP <- score_summary$n[1]

top_scoring <- samples@meta.data %>%
  slice_max(order_by = PolyKRT_module_score, n = N_TOP)

concentration <- top_scoring %>%
  count(cell_type_unintegrated, sort = TRUE) %>%
  mutate(pct_of_top_group = round(100 * n / N_TOP, 1))

cat("\n=== Cell type composition of the top", N_TOP,
    "highest PolyKRT-scoring cells (unintegrated) ===\n")
print(concentration)

write.csv(concentration,
          paste0(DATE_PREFIX, "-PolyKRT-module-concentration-test-unintegrated.csv"),
          row.names = FALSE)

cat("\nINTERPRETATION:\n")
cat("  If top scorers concentrate heavily (e.g. >70%) in ONE cell type, the\n")
cat("  signature maps to a distinct, identifiable population in this data.\n")
cat("  If spread broadly across multiple types, it does not cleanly demarcate\n")
cat("  a single population here -- report honestly either way.\n")





# -----------------------------------------------------------------------------
# Bautista immature TEC marker heatmap -- unintegrated, subset
# (Bipolar TEC, Pericyte, Fibroblast, Erythroid excluded)
# -----------------------------------------------------------------------------
exclude_types <- c("Bipolar TEC", "Pericyte", "Fibroblast", "Erythroid")

samples_vln_subset <- subset(
  samples,
  subset = !(cell_type_unintegrated %in% exclude_types)
)
samples_vln_subset$cell_type_unintegrated <- droplevels(samples_vln_subset$cell_type_unintegrated)

cat("Cell types remaining:\n")
print(table(samples_vln_subset$cell_type_unintegrated))
stopifnot(!any(levels(samples_vln_subset$cell_type_unintegrated) %in% exclude_types))

cluster_order_subset <- c("Nurse", "cTEC", "BMC TEC", "mTEC lo", "mTEC hi",
                          "Mimetic", "Endothelial")

stopifnot(setequal(cluster_order_subset,
                   levels(samples_vln_subset$cell_type_unintegrated)))

avg_mat_subset <- AverageExpression(
  samples_vln_subset, features = bautista_genes, assay = "RNA",
  group.by = "cell_type_unintegrated", verbose = TRUE
)$RNA
avg_mat_subset <- avg_mat_subset[bautista_genes, cluster_order_subset, drop = FALSE]

log_mat_subset <- log1p(avg_mat_subset)
z_mat_subset   <- t(scale(t(log_mat_subset)))
z_mat_subset[is.nan(z_mat_subset)] <- 0

# Sanity check
stopifnot(nrow(z_mat_subset) == length(bautista_genes))
stopifnot(ncol(z_mat_subset) == length(cluster_order_subset))
stopifnot(!any(is.na(z_mat_subset)))

df_heatmap_subset <- as.data.frame(z_mat_subset) %>%
  tibble::rownames_to_column("gene") %>%
  tidyr::pivot_longer(-gene, names_to = "cluster", values_to = "z") %>%
  mutate(
    gene    = factor(gene, levels = rev(bautista_genes)),
    cluster = factor(cluster, levels = cluster_order_subset)
  )

p_bautista_heatmap_subset <- ggplot(df_heatmap_subset, aes(x = cluster, y = gene, fill = z)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient2(
    low = "#F7E27A", mid = "grey95", high = "#1B3B6F",
    midpoint = 0, name = "Relative\nexpr. (z)"
  ) +
  labs(title = "Select Markers") +
  theme_minimal(base_size = 10) +
  theme(
    axis.title  = element_blank(),
    plot.title  = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid  = element_blank()
  )
print(p_bautista_heatmap_subset)
ggsave(paste0(DATE_PREFIX, "-bautista-immature-TEC-heatmap-unintegrated-subset.pdf"),
       p_bautista_heatmap_subset, width = 5, height = 4)





# -----------------------------------------------------------------------------
# Nurse subcluster characterization -- clusters 3, 13, 15 (unintegrated)
# Marker enrichment (already computed in markers_unint) + cell cycle scoring
# -----------------------------------------------------------------------------
markers_unint <- read.csv(paste0("results/260805-2000-unintegrated-markers-FULL.csv"))

top15_unint <- markers_unint %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 15) %>%
  ungroup()

# Sanity check
stopifnot(nrow(markers_unint) > 0)
cat("Clusters with markers:", length(unique(markers_unint$cluster)), "\n")

# nurse cluster numbers in unintegrated data
nurse_clusters <- c("3", "13", "15")

# Cell cycle scoring, if not already run on this object
samples <- CellCycleScoring(
  samples,
  s.features   = cc.genes.updated.2019$s.genes,
  g2m.features = cc.genes.updated.2019$g2m.genes,
  set.ident    = FALSE, verbose = TRUE
)
stopifnot(!any(is.na(samples$Phase)))





# -----------------------------------------------------------------------------
# Marker genes per Nurse subcluster (already computed, just pull and label)
# -----------------------------------------------------------------------------
nurse_markers_top15 <- top15_unint %>%
  filter(cluster %in% nurse_clusters) %>%
  select(cluster, gene, avg_log2FC, pct.1, pct.2, p_val_adj)

cat("=== Top markers per Nurse subcluster ===\n")
print(nurse_markers_top15, n = 45)

write.csv(nurse_markers_top15,
          paste0(DATE_PREFIX, "-nurse-subcluster-markers-unintegrated.csv"),
          row.names = FALSE)





# -----------------------------------------------------------------------------
# Cell cycle phase per Nurse subcluster
# -----------------------------------------------------------------------------
nurse_cellcycle <- samples@meta.data %>%
  filter(clusters_unintegrated %in% nurse_clusters) %>%
  group_by(clusters_unintegrated, Phase) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Phase, values_from = n, values_fill = 0) %>%
  mutate(
    total = G1 + S + G2M,
    pct_cycling = round(100 * (S + G2M) / total, 1)
  )

cat("\n=== Cell cycle phase by Nurse subcluster ===\n")
print(nurse_cellcycle)

write.csv(nurse_cellcycle,
          paste0(DATE_PREFIX, "-nurse-subcluster-cellcycle-unintegrated.csv"),
          row.names = FALSE)

# Sanity check
stopifnot(sum(nurse_cellcycle$total) == sum(samples$clusters_unintegrated %in% nurse_clusters))





# -----------------------------------------------------------------------------
# Heatmap: top 15 markers defining Bipolar TEC (unintegrated)
# -----------------------------------------------------------------------------
bipolar_genes <- top15_unint %>%
  filter(cluster == "7") %>%
  arrange(desc(avg_log2FC)) %>%
  pull(gene)

stopifnot(length(bipolar_genes) == 15)
cat("Bipolar TEC top 15 genes:\n")
cat(paste(bipolar_genes, collapse = ", "), "\n")

exclude_types <- c("Pericyte", "Fibroblast", "Erythroid")   # Bipolar TEC retained

samples_bipolar_subset <- subset(
  samples,
  subset = !(cell_type_unintegrated %in% exclude_types)
)
samples_bipolar_subset$cell_type_unintegrated <- droplevels(samples_bipolar_subset$cell_type_unintegrated)

cat("\nCell types remaining:\n")
print(table(samples_bipolar_subset$cell_type_unintegrated))
stopifnot(!any(levels(samples_bipolar_subset$cell_type_unintegrated) %in% exclude_types))

cluster_order_bipolar <- c("Nurse", "cTEC", "BMC TEC", "mTEC lo", "mTEC hi",
                           "Mimetic","Endothelial", "Bipolar TEC")

stopifnot(setequal(cluster_order_bipolar,
                   levels(samples_bipolar_subset$cell_type_unintegrated)))

avg_mat_bipolar <- AverageExpression(
  samples_bipolar_subset, features = bipolar_genes, assay = "RNA",
  group.by = "cell_type_unintegrated", verbose = TRUE
)$RNA
avg_mat_bipolar <- avg_mat_bipolar[bipolar_genes, cluster_order_bipolar, drop = FALSE]

log_mat_bipolar <- log1p(avg_mat_bipolar)
z_mat_bipolar   <- t(scale(t(log_mat_bipolar)))
z_mat_bipolar[is.nan(z_mat_bipolar)] <- 0

# Sanity check
stopifnot(nrow(z_mat_bipolar) == 15)
stopifnot(ncol(z_mat_bipolar) == length(cluster_order_bipolar))
stopifnot(!any(is.na(z_mat_bipolar)))

df_heatmap_bipolar <- as.data.frame(z_mat_bipolar) %>%
  tibble::rownames_to_column("gene") %>%
  tidyr::pivot_longer(-gene, names_to = "cluster", values_to = "z") %>%
  mutate(
    gene    = factor(gene, levels = rev(bipolar_genes)),
    cluster = factor(cluster, levels = cluster_order_bipolar)
  )

p_bipolar_heatmap <- ggplot(df_heatmap_bipolar, aes(x = cluster, y = gene, fill = z)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient2(
    low = "#F7E27A", mid = "grey95", high = "#1B3B6F",
    midpoint = 0, name = "Relative\nexpr. (z)"
  ) +
  labs(title = "Select Markers") +
  theme_minimal(base_size = 10) +
  theme(
    axis.title  = element_blank(),
    plot.title  = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid  = element_blank()
  )
print(p_bipolar_heatmap)
ggsave(paste0(DATE_PREFIX, "-bipolar-TEC-marker-heatmap-unintegrated.pdf"),
       p_bipolar_heatmap, width = 7, height = 6)





# -----------------------------------------------------------------------------
# Bipolar TEC subcluster analysis -- Step 1: visualize in isolation
# -----------------------------------------------------------------------------
bipolar_cells <- subset(samples, subset = clusters_unintegrated == "7")

cat("Bipolar TEC total cells:", ncol(bipolar_cells), "\n")

p_bipolar_alone <- DimPlot(bipolar_cells, reduction = "umap_unintegrated",
                           group.by = "clusters_unintegrated") +
  NoLegend() +
  labs(title = "Bipolar TEC only (cluster 7, unintegrated)")
print(p_bipolar_alone)
ggsave(paste0(DATE_PREFIX, "-bipolar-TEC-alone-umap.pdf"), p_bipolar_alone,
       width = 6, height = 5)

# Also on the full UMAP, highlighted against everything else -- useful for
# seeing exactly which neighbors it's spatially adjacent to
p_bipolar_highlight <- DimPlot(
  samples, reduction = "umap_unintegrated",
  cells.highlight = WhichCells(samples, expression = clusters_unintegrated == "7")
) +
  labs(title = "Bipolar TEC highlighted on full UMAP") +
  scale_color_manual(values = c("grey80", "red"), labels = c("Other", "Bipolar TEC"))
print(p_bipolar_highlight)
ggsave(paste0(DATE_PREFIX, "-bipolar-TEC-highlighted-full-umap.pdf"),
       p_bipolar_highlight, width = 8, height = 7)


c7_coords <- Embeddings(samples, "umap_unintegrated") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell") %>%
  rename(umap_1 = 2, umap_2 = 3)
c7_coords$clusters_unintegrated <- as.character(samples$clusters_unintegrated[c7_coords$cell])
c7_coords <- c7_coords %>% filter(clusters_unintegrated == "7")

c7_coords <- c7_coords %>%
  mutate(sub_region = case_when(
    umap_1 >= -7 & umap_1 <= 0  & umap_2 >= 0  & umap_2 <= 6   ~ "near_cTEC",
    umap_1 >= 5  & umap_1 <= 9  & umap_2 >= 5  & umap_2 <= 8.5 ~ "near_Nurse",
    umap_1 >= -6 & umap_1 <= -1 & umap_2 >= -9 & umap_2 <= -4  ~ "near_mTEC_lo",
    umap_1 >= -1 & umap_1 <= 5  & umap_2 >= -9 & umap_2 <= -4  ~ "near_mTEC_hi",
    TRUE ~ "unassigned"
  ))

cat("=== Bipolar TEC sub-region cell counts ===\n")
print(table(c7_coords$sub_region))
stopifnot(sum(table(c7_coords$sub_region)) == nrow(c7_coords))


# Attach sub-region labels back to the full samples object
samples$c7_subregion <- NA_character_
samples$c7_subregion[match(c7_coords$cell, colnames(samples))] <- c7_coords$sub_region

# Sanity check
stopifnot(sum(!is.na(samples$c7_subregion)) == nrow(c7_coords))

# Plot: Bipolar TEC only, colored by sub-region assignment
bipolar_cells <- subset(samples, subset = clusters_unintegrated == "7")

p_c7_gating <- DimPlot(bipolar_cells, reduction = "umap_unintegrated",
                       group.by = "c7_subregion", label = TRUE, repel = TRUE) +
  labs(title = "Bipolar TEC gated into 4 sub-regions (unintegrated)")
print(p_c7_gating)
ggsave(paste0(DATE_PREFIX, "-bipolar-TEC-subregion-gating.pdf"), p_c7_gating,
       width = 7, height = 6)





# -----------------------------------------------------------------------------
# Bipolar TEC sub-region marker analysis: one-vs-all + one-vs-specific-neighbor
# -----------------------------------------------------------------------------
neighbor_map <- list(
  near_cTEC    = c("0", "2"),
  near_mTEC_lo = c("1"),
  near_mTEC_hi = c("6"),
  near_Nurse   = c("3", "13", "15")
)

# Sanity check: confirm every listed neighbor cluster exists and has cells
for (region in names(neighbor_map)) {
  ids <- neighbor_map[[region]]
  missing_ids <- setdiff(ids, unique(as.character(samples$clusters_unintegrated)))
  if (length(missing_ids) > 0) {
    stop(region, ": cluster ID(s) not found: ", paste(missing_ids, collapse = ", "))
  }
  n_neighbor_cells <- sum(samples$clusters_unintegrated %in% ids)
  cat(sprintf("%-14s neighbor clusters %-10s -> %d cells\n",
              region, paste(ids, collapse = ","), n_neighbor_cells))
}

# Give every cell a per-subregion identity for marker testing
samples$ident_for_c7_analysis <- ifelse(
  !is.na(samples$c7_subregion) & samples$c7_subregion != "unassigned",
  samples$c7_subregion,
  as.character(samples$clusters_unintegrated)
)
Idents(samples) <- "ident_for_c7_analysis"
stopifnot(identical(as.character(Idents(samples)), samples$ident_for_c7_analysis))

results_all      <- list()
results_neighbor <- list()

for (region in names(neighbor_map)) {
  cat("\n", strrep("=", 60), "\n")
  cat("SUB-REGION:", region, "\n")
  cat(strrep("=", 60), "\n")
  
  n_cells <- sum(samples$ident_for_c7_analysis == region, na.rm = TRUE)
  cat("Cell count:", n_cells, "\n")
  if (n_cells < 20) cat("WARNING: few cells -- results will be noisy.\n")
  
  m_all <- tryCatch({
    FindMarkers(samples, ident.1 = region, only.pos = TRUE,
                min.pct = 0.25, logfc.threshold = 0.25, verbose = TRUE)
  }, error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); NULL })
  
  if (!is.null(m_all)) {
    m_all$gene <- rownames(m_all)
    m_all <- m_all %>% arrange(desc(avg_log2FC))
    results_all[[region]] <- m_all
    cat("\nTop 15, vs-ALL:\n")
    print(head(m_all[, c("gene","avg_log2FC","pct.1","pct.2","p_val_adj")], 15))
    write.csv(m_all, paste0(DATE_PREFIX, "-c7-", region, "-vs-ALL.csv"), row.names = FALSE)
  }
  
  m_neighbor <- tryCatch({
    FindMarkers(samples, ident.1 = region, ident.2 = neighbor_map[[region]],
                only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = TRUE)
  }, error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); NULL })
  
  if (!is.null(m_neighbor)) {
    m_neighbor$gene <- rownames(m_neighbor)
    m_neighbor <- m_neighbor %>% arrange(desc(avg_log2FC))
    results_neighbor[[region]] <- m_neighbor
    cat("\nTop 15, vs-SPECIFIC-NEIGHBOR (", paste(neighbor_map[[region]], collapse=","), "):\n")
    print(head(m_neighbor[, c("gene","avg_log2FC","pct.1","pct.2","p_val_adj")], 15))
    write.csv(m_neighbor, paste0(DATE_PREFIX, "-c7-", region, "-vs-NEIGHBOR.csv"), row.names = FALSE)
  }
}





# -----------------------------------------------------------------------------
# Sanity check: confirm c7 sub-region analysis used unintegrated cluster IDs
# only, not harmony cluster IDs, at every step (gating, neighbor_map, Idents)
# -----------------------------------------------------------------------------

# 1. Bipolar TEC cluster 7 was gated from clusters_unintegrated, not clusters_harmony
stopifnot(all(samples$clusters_unintegrated[!is.na(samples$c7_subregion)] == "7"))
cat("c7_subregion cells are ALL clusters_unintegrated == 7: TRUE\n")

# 2. neighbor_map cluster IDs must exist in clusters_unintegrated, and be
#    checked against clusters_harmony to confirm they're not silently valid
#    there too by coincidence (which would mask a branch mix-up)
for (region in names(neighbor_map)) {
  ids <- neighbor_map[[region]]
  in_unint   <- all(ids %in% levels(samples$clusters_unintegrated))
  cat(sprintf("%-14s ids=%-10s in unintegrated levels: %s\n",
              region, paste(ids, collapse=","), in_unint))
}
stopifnot(all(sapply(neighbor_map, function(ids) all(ids %in% levels(samples$clusters_unintegrated)))))

# 3. Confirm Idents used for FindMarkers were built from the unintegrated-only
#    ident_for_c7_analysis column, not clusters_harmony, immediately before the loop
stopifnot(identical(as.character(Idents(samples)), samples$ident_for_c7_analysis))
stopifnot(!identical(samples$ident_for_c7_analysis, as.character(samples$clusters_harmony)))
cat("Idents at time of FindMarkers matched ident_for_c7_analysis (unintegrated-derived): TRUE\n")

# 4. Cross-tab: how do the 4 sub-regions distribute across clusters_harmony?
#    (informative only -- NOT expected to equal neighbor_map, since harmony
#    has different cluster numbering/composition; just confirms no accidental
#    harmony-cluster-number reuse crept into the neighbor_map ids)
print(table(samples$c7_subregion, samples$clusters_harmony, useNA = "no"))






# -----------------------------------------------------------------------------
# Bipolar TEC sub-region analysis, Approach B: internal distinctiveness,
# then cross-reference against actual neighbor cluster markers
# -----------------------------------------------------------------------------

# Explicitly set Idents -- do not assume prior state
Idents(samples) <- "ident_for_c7_analysis"
stopifnot(all(as.character(Idents(samples)) == unname(samples$ident_for_c7_analysis)))

subregions <- c("near_cTEC", "near_mTEC_lo", "near_mTEC_hi", "near_Nurse")

results_internal <- list()

for (region in subregions) {
  stopifnot(all(as.character(Idents(samples)) == unname(samples$ident_for_c7_analysis)))
  
  other_regions <- setdiff(subregions, region)
  cat("\n", strrep("=", 60), "\n")
  cat("SUB-REGION:", region, "vs. other 3 sub-regions (internal)\n")
  cat(strrep("=", 60), "\n")
  
  n_cells <- sum(samples$ident_for_c7_analysis == region, na.rm = TRUE)
  cat("Cell count:", n_cells, "\n")
  
  m <- tryCatch({
    FindMarkers(samples, ident.1 = region, ident.2 = other_regions,
                only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,
                verbose = TRUE)
  }, error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); NULL })
  
  if (!is.null(m)) {
    m$gene <- rownames(m)
    m <- m %>% arrange(desc(avg_log2FC))
    results_internal[[region]] <- m
    cat("\nTop 15, internal distinctiveness:\n")
    print(head(m[, c("gene","avg_log2FC","pct.1","pct.2","p_val_adj")], 15))
    write.csv(m, paste0(DATE_PREFIX, "-c7-", region, "-vs-OTHER3-internal.csv"),
              row.names = FALSE)
  }
}

# -----------------------------------------------------------------------------
# Cross-reference: does each sub-region's internal marker list overlap with
# its ACTUAL neighbor's own top markers?
# -----------------------------------------------------------------------------
neighbor_top_markers <- list(
  near_cTEC    = top15_unint %>% filter(cluster %in% c("0","2")) %>% pull(gene) %>% unique(),
  near_mTEC_lo = top15_unint %>% filter(cluster == "1") %>% pull(gene) %>% unique(),
  near_mTEC_hi = top15_unint %>% filter(cluster == "6") %>% pull(gene) %>% unique(),
  near_Nurse   = top15_unint %>% filter(cluster %in% c("3","13","15")) %>% pull(gene) %>% unique()
)

cat("\n", strrep("=", 60), "\n")
cat("CROSS-REFERENCE: internal markers vs. actual neighbor's top markers\n")
cat(strrep("=", 60), "\n")

for (region in subregions) {
  if (is.null(results_internal[[region]])) next
  
  internal_top20 <- results_internal[[region]] %>%
    slice_max(order_by = avg_log2FC, n = 20) %>%
    pull(gene)
  
  neighbor_genes <- neighbor_top_markers[[region]]
  overlap <- intersect(internal_top20, neighbor_genes)
  
  cat("\n", region, ":\n", sep = "")
  cat("  Internal top-20 genes:", paste(internal_top20, collapse = ", "), "\n")
  cat("  Neighbor's own top genes:", paste(neighbor_genes, collapse = ", "), "\n")
  cat("  OVERLAP (", length(overlap), " genes):", paste(overlap, collapse = ", "), "\n", sep = "")
}





# =============================================================================
# Bipolar TEC standalone re-embedding: subset, wipe stale reductions,
# fresh PCA/UMAP/clustering on just these cells. c7_subregion carried along
# as a cross-reference label only, NOT used to build the embedding.
# =============================================================================

# -----------------------------------------------------------------------------
# Step 1: Subset to Bipolar TEC cells only
# -----------------------------------------------------------------------------
bipolar_standalone <- subset(samples, subset = clusters_unintegrated == "7")

cat("Bipolar TEC standalone cell count:", ncol(bipolar_standalone), "\n")
stopifnot(ncol(bipolar_standalone) == 1264)

cat("\nc7_subregion distribution (carried along, not used for clustering):\n")
print(table(bipolar_standalone$c7_subregion, useNA = "ifany"))

# -----------------------------------------------------------------------------
# Step 2: Wipe stale layers/reductions from the parent object
# -----------------------------------------------------------------------------
DefaultAssay(bipolar_standalone) <- "RNA"
bipolar_standalone <- DietSeurat(
  bipolar_standalone,
  layers    = c("counts", "data"),
  dimreducs = NULL,
  graphs    = NULL,
  misc      = FALSE
)

# Sanity check: confirm the wipe worked, c7_subregion metadata survived
stopifnot(length(Reductions(bipolar_standalone)) == 0)
stopifnot("c7_subregion" %in% colnames(bipolar_standalone@meta.data))
cat("\nLayers after wipe:", paste(Layers(bipolar_standalone, assay = "RNA"), collapse = ", "), "\n")

# -----------------------------------------------------------------------------
# Step 3: Fresh preprocessing on just these 1,264 cells
# -----------------------------------------------------------------------------
bipolar_standalone <- NormalizeData(bipolar_standalone, verbose = TRUE)
bipolar_standalone <- FindVariableFeatures(bipolar_standalone, nfeatures = 2000, verbose = TRUE)
bipolar_standalone <- ScaleData(bipolar_standalone, verbose = TRUE)
bipolar_standalone <- RunPCA(bipolar_standalone, npcs = 30, verbose = TRUE)

p_elbow <- ElbowPlot(bipolar_standalone, ndims = 30) +
  labs(title = "Bipolar TEC standalone PCA: elbow plot")
print(p_elbow)
ggsave(paste0(DATE_PREFIX, "-bipolar-standalone-elbow.pdf"), p_elbow, width = 6, height = 5)

# -----------------------------------------------------------------------------
# Step 4: Neighbors, clustering (sweep a few resolutions), UMAP
# -----------------------------------------------------------------------------
N_PCS <- 10

bipolar_standalone <- FindNeighbors(bipolar_standalone, dims = 1:N_PCS, verbose = TRUE)

resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0)
for (r in resolutions) {
  bipolar_standalone <- FindClusters(bipolar_standalone, resolution = r, verbose = TRUE)
}

bipolar_standalone <- RunUMAP(bipolar_standalone, dims = 1:N_PCS, verbose = TRUE)

cat("\n=== Cluster counts at each resolution ===\n")
for (r in resolutions) {
  col <- paste0("RNA_snn_res.", r)
  cat(sprintf("res %.1f -> %d clusters\n", r, length(unique(bipolar_standalone[[col]][[1]]))))
}

# check number of cells in each cluster across the resolution sweep
for (r in resolutions) {
  col <- paste0("RNA_snn_res.", r)
  cat("\n--- res =", r, "---\n")
  print(table(bipolar_standalone[[col]][[1]]))
  }

# -----------------------------------------------------------------------------
# Step 5: Visualize resolutions, colored by data-driven cluster
# -----------------------------------------------------------------------------
p_res <- wrap_plots(
  lapply(resolutions, function(r) {
    col <- paste0("RNA_snn_res.", r)
    DimPlot(bipolar_standalone, group.by = col, label = TRUE) +
      labs(title = paste0("res = ", r)) + NoLegend()
  }),
  ncol = 3
)
print(p_res)
ggsave(paste0(DATE_PREFIX, "-bipolar-standalone-resolutions.pdf"), p_res,
       width = 12, height = 8)

# -----------------------------------------------------------------------------
# Step 6: Cross-reference -- does the DATA-DRIVEN structure line up with the
# original UMAP-coordinate-gated c7_subregion labels?
# -----------------------------------------------------------------------------
# Do this at a couple of resolutions once you've picked one from Step 5 --
# placeholder shown at res=0.4, adjust after visual inspection
CHOSEN_RES <- 0.2
res_col <- paste0("RNA_snn_res.", CHOSEN_RES)

cat("\n=== Data-driven subcluster (res=", CHOSEN_RES,
    ") vs. coordinate-gated c7_subregion ===\n", sep = "")
cross_tab <- table(bipolar_standalone[[res_col]][[1]], bipolar_standalone$c7_subregion)
print(cross_tab)

# Sanity check
stopifnot(sum(cross_tab) == ncol(bipolar_standalone) - sum(is.na(bipolar_standalone$c7_subregion)))

# Side-by-side visual: data-driven clusters vs. coordinate-gated regions,
# same UMAP embedding
p_compare <- (
  DimPlot(bipolar_standalone, group.by = res_col, label = TRUE) +
    labs(title = paste0("Data-driven clusters (res=", CHOSEN_RES, ")"))
) | (
  DimPlot(bipolar_standalone, group.by = "c7_subregion", label = TRUE) +
    labs(title = "Coordinate-gated sub-regions (for comparison)")
)
print(p_compare)
ggsave(paste0(DATE_PREFIX, "-bipolar-standalone-vs-coordinate-gating.pdf"),
       p_compare, width = 12, height = 6)

saveRDS(bipolar_standalone, paste0(DATE_PREFIX, "-bipolar-TEC-standalone.rds"))

























# highlight individual clusters overlay over the entire umap
highlight_cluster <- function(cluster_num) {
  DimPlot(samples, reduction = "umap_unintegrated",
          cells.highlight = WhichCells(samples, expression = clusters_unintegrated == as.character(cluster_num))) +
    scale_color_manual(values = c("grey80", "red"), labels = c("Other", paste0("Cluster ", cluster_num))) +
    labs(title = paste0("Cluster ", cluster_num, " highlighted"))
}

highlight_cluster(0)
highlight_cluster(1)
highlight_cluster(2)
highlight_cluster(3)
highlight_cluster(4)
highlight_cluster(5)
highlight_cluster(6)
highlight_cluster(7)
highlight_cluster(8)
highlight_cluster(9)
highlight_cluster(10)
highlight_cluster(11)
highlight_cluster(12)
highlight_cluster(13)
highlight_cluster(14)
highlight_cluster(15)
highlight_cluster(16)
highlight_cluster(17)
highlight_cluster(18)
highlight_cluster(19)
highlight_cluster(20)
highlight_cluster(21)
highlight_cluster(22)
highlight_cluster(23)
highlight_cluster(24)
highlight_cluster(25)

# visualize clusters on umap 
DimPlot(samples, reduction = "umap_unintegrated",
        group.by = "clusters_unintegrated", label = TRUE) +
  labs(title = "Unintegrated: clusters")

DimPlot(samples, reduction = "umap_harmony",
        group.by = "clusters_harmony", label = TRUE) +
  labs(title = "Harmony-integrated: clusters")

DimPlot(
  samples, reduction = "umap_unintegrated", label = TRUE,
  group.by = "clusters_unintegrated",
  cells.highlight = WhichCells(samples, expression = clusters_unintegrated == "7")
) +
  labs(title = "Unintegrated: cluster 7 highlighted") +
  scale_color_manual(values = c("grey80", "red"), labels = c("Other", "Cluster 7"))

DimPlot(
  samples, reduction = "umap_harmony",
  cells.highlight = WhichCells(samples, expression = clusters_harmony == "6")
) +
  labs(title = "Harmony: cluster 6 highlighted") +
  scale_color_manual(values = c("grey80", "red"), labels = c("Other", "Cluster 6"))

p_full_reference <- DimPlot(samples, reduction = "umap_unintegrated",
                            group.by = "cell_type_unintegrated", label = TRUE, repel = TRUE) +
  labs(title = "Full annotated UMAP (unintegrated) -- for spatial reference")
print(p_full_reference)


# -----------------------------------------------------------------------------
# Doublet detection and removal -- PER SAMPLE, before merge
# -----------------------------------------------------------------------------
remove_doublets <- function(seurat_obj, sample_name) {
  cat("\nRunning scDblFinder on:", sample_name, "\n")
  sce <- as.SingleCellExperiment(seurat_obj, assay = "RNA")
  sce <- scDblFinder(sce, verbose = TRUE)
  
  seurat_obj$scDblFinder_class <- as.character(sce$scDblFinder.class)
  seurat_obj$scDblFinder_score <- as.numeric(sce$scDblFinder.score)
  
  n_before <- ncol(seurat_obj)
  seurat_obj <- subset(seurat_obj, subset = scDblFinder_class == "singlet")
  n_after <- ncol(seurat_obj)
  
  cat(sample_name, ": ", n_before, " -> ", n_after,
      " (", round(100 * (n_before - n_after) / n_before, 1), "% doublets removed)\n", sep = "")
  seurat_obj
}

ht67_cd205neg_sub <- remove_doublets(ht67_cd205neg_sub, "ht67-cd205neg")
ht67_cd205pos_sub <- remove_doublets(ht67_cd205pos_sub, "ht67-cd205pos")
ht70_sub           <- remove_doublets(ht70_sub, "ht70")
ht71_sub           <- remove_doublets(ht71_sub, "ht71")

sanity_cells <- c(ht67_cd205neg_sub = ncol(ht67_cd205neg_sub),
                  ht67_cd205pos_sub = ncol(ht67_cd205pos_sub),
                  ht70_sub          = ncol(ht70_sub),
                  ht71_sub          = ncol(ht71_sub))
cat("\nCells remaining per sample after doublet removal:\n")
print(sanity_cells)
stopifnot(all(sanity_cells > 0))

