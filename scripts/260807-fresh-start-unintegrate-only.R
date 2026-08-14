# 10x flex data set of human thymic epithellial cells from samples HT67, HT70, HT71

library(Seurat)
library(harmony)
library(scDblFinder)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(patchwork)
library(scCustomize)
library(ggrastr)

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
# Find Neighbor, Clusters, UMAP
# -----------------------------------------------------------------------------
samples <- FindNeighbors(samples, reduction = "pca", dims = 1:18,
                         graph.name = c("nn_unintegrated", "snn_unintegrated"))
samples <- FindClusters(samples, graph.name = "snn_unintegrated",
                        resolution = 0.5, cluster.name = "clusters_unintegrated")
samples <- RunUMAP(samples, reduction = "pca", dims = 1:18,
                   reduction.name = "umap_unintegrated")

# -----------------------------------------------------------------------------
# Sanity checks
# -----------------------------------------------------------------------------
cat("Reductions present:", paste(Reductions(samples), collapse = ", "), "\n")
cat("Unintegrated clusters found:", length(unique(samples$clusters_unintegrated)), "\n")
stopifnot(all(c("umap_unintegrated") %in% Reductions(samples)))
stopifnot(!any(is.na(samples$clusters_unintegrated)))

table(samples$clusters_unintegrated)

#saveRDS(samples, paste0(DATE_PREFIX, "-samples-both-pre-and-post-integration.rds"))
#samples <- readRDS(file = "data/rds-objects/260807-1222-samples-both-pre-and-post-integration.rds")

# -----------------------------------------------------------------------------
# Find markers
# -----------------------------------------------------------------------------
Idents(samples) <- "clusters_unintegrated"
stopifnot(identical(as.character(Idents(samples)), as.character(samples$clusters_unintegrated)))
markers_unint <- FindAllMarkers(samples, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Full marker tables
write.csv(markers_unint, paste0(DATE_PREFIX, "-unintegrated-markers-FULL.csv"), row.names = FALSE)

# Top 15 per cluster, for a quick scan
top15_unint <- markers_unint %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 15) %>%
  ungroup()

write.csv(top15_unint, paste0(DATE_PREFIX, "-unintegrated-markers-TOP15.csv"), row.names = FALSE)

saveRDS(samples, paste0(DATE_PREFIX, "-fresh-start-unintegrate-only-post-find-markers.rds"))
#samples <- readRDS(file = "data/rds-objects/260807-1222-fresh-start-unintegrate-only-post-find-markers.rds")

# -----------------------------------------------------------------------------
# Annotating cell types
# -----------------------------------------------------------------------------
set.seed(2)

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
  "cTEC",               # this cTEC cluster contains some mesenchymal signatures    
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

# Final sanity check
stopifnot(!any(is.na(samples$cell_type_unintegrated)))
cat("\n=== Final annotation counts ===\n")
print(table(samples$cell_type_unintegrated))

# check annotation on umap
DimPlot(samples, reduction = "umap_unintegrated",
        group.by = "cell_type_unintegrated", label = TRUE, repel = TRUE) +
  labs(title = "Unintegrated: annotated cell types")

DimPlot(samples, reduction = "umap_unintegrated", 
        group.by = "cell_type_unintegrated",
        label = FALSE,
        cols = c("#D8A767","#F47D2B","#D24B27","#41B6C4","#E7298A",
                 "#0570B0","#3D3D3D","#89C75F", "#7E1416" , "#E6C2DC"))

saveRDS(samples, paste0(DATE_PREFIX, "-samples-unintegrated-only-post-annotating.rds"))
#samples <- readRDS(file = "data/rds-objects/260812-1100-samples-unintegrated-only-post-annotating.rds")

# -----------------------------------------------------------------------------
# DimPlot
# -----------------------------------------------------------------------------
exclude_types_dimplot <- c("Bipolar TEC","Endothelial","Pericyte","Erythroid")

samples_dimplot_subset <- subset(
  samples,
  subset = !(cell_type_unintegrated %in% exclude_types_dimplot)
)
samples_dimplot_subset$cell_type_unintegrated <- droplevels(samples_dimplot_subset$cell_type_unintegrated)

# Sanity check
cat("Cell types remaining:\n")
print(table(samples_dimplot_subset$cell_type_unintegrated))
stopifnot(!any(levels(samples_dimplot_subset$cell_type_unintegrated) %in% exclude_types_dimplot))

p_dimplot_subset <- DimPlot(
  samples_dimplot_subset, reduction = "umap_unintegrated",
  group.by = "cell_type_unintegrated", label = FALSE, repel = TRUE,
  cols = c("#D8A767","#208A42","#D24B27","#41B6C4","#E7298A","#0570B0",
           "#89C75F"),
  split.by = "sample.id"
) +
  labs(title = "Unintegrated: annotated cell types (subset)")
print(p_dimplot_subset)
ggsave(paste0(DATE_PREFIX, "-dimplot-unintegrated-subset.pdf"),
       p_dimplot_subset, width = 24, height = 7)

palman_values =c("#D8A767","#3D3D3D","#D51F26","#F47D2B","#0C727C","#3BBCA8",
                 "#FEE500","#D24B27","#89288F","#9983BD","#7E1416","#272E6A",
                 "#F37B7D","#89C75F","#208A42","#90D5E4","#E6C2DC","#E7298A",
                 "#0570B0","#A8DDB5","#FEB24C","#41B6C4","#BCBDDC")


# batch save featureplots of interest
save_features <- function(feature, seur, reduction = "umap-unintegrated", label = TRUE, 
                          label_size = 4, order = TRUE, split.by = NULL) {
  
  p <- FeaturePlot_scCustom(seur, reduction = reduction, features = feature,
                            colors_use = viridis_light_high, label = label, 
                            label.size = label_size, order = order, 
                            split.by = split.by)
  
  ggsave(plot=p, paste0("results/",DATE_PREFIX,"-featureplot-",obj_name,feature, ".pdf"),
         width=8, height=7)
}

foi_nurse <- c("TBATA","PRSS16","PSMB11","PTPRC","CD3E")
foi_ctec <- c("TBATA","PRSS16","PSMB11","LY75")
foi_bmc_tec <- c("COL4A5","COL4A6","COL8A1","ITGB4")
foi_mtec_lo <- c("CCL19","KRT15","KRT19","FN1","CLU")
foi_mtec_hi <- c("AIRE","FEZF2")
foi_mimetic <- "CD24"
foi_endothelial <- "PECAM1"
foi_all_tec <- "TP63"
foi_bautista_immature_tec_1 <- c("ACTB","JUNB","FOS")
foi_bautista_immature_tec_2 <- c("IGFBP5","NNMT","MAOA","DPYS","FKBP5","GLUL")
foi_bautista_mtec_lo <- c("GABRA5","LYPD1")
foi_other_markers <- c("CD200","ITGA6","BCAM","EPCAM","CD24","LY75","ITGB4",
                       "CD74")

save_features <- function(feature, seur, reduction = "umap_unintegrated",
                          label = FALSE, label_size = 4, 
                          order = TRUE, split.by = NULL) {
  
  p <- FeaturePlot_scCustom(seur, reduction = reduction, features = feature,
                            colors_use = viridis_light_high, label = label, 
                            label.size = label_size, order = order, 
                            split.by = split.by)
  
  ggsave(plot=p, paste0("results/",DATE_PREFIX,"-featureplot-",feature, ".pdf"),
         width=8, height=7)
}

sapply(foi_nurse, save_features,  seur=samples)
sapply(foi_ctec, save_features,  seur=samples)
sapply(foi_bmc_tec, save_features,  seur=samples)
sapply(foi_mtec_lo, save_features,  seur=samples)
sapply(foi_mtec_hi, save_features,  seur=samples)
sapply(foi_mimetic, save_features,  seur=samples)
sapply(foi_endothelial, save_features,  seur=samples)
sapply(foi_all_tec, save_features,  seur=samples)
sapply(foi_bautista_immature_tec_1, save_features,  seur=samples)
sapply(foi_bautista_immature_tec_2, save_features,  seur=samples)
sapply(foi_bautista_mtec_lo, save_features,  seur=samples)
sapply(foi_other_markers, save_features,  seur=samples)

# -----------------------------------------------------------------------------
# "Select Markers" heatmap -- unintegrated, 7-category subset
# (Bipolar TEC, Pericyte, Erythroid excluded
# -----------------------------------------------------------------------------
select_markers_genes <- c(
  "PTPRC", "PDPN",
  "LY75", "PSMB11", "PRSS16",
  "COL4A5", "COL4A6","COL8A1","ITGB4",
  "CD200", "CLU","AIRE", "FEZF2",
  "CD24","ITGA6", "BCAM","CD74"
  # "PECAM1", "PDGFRB","RGS5",
  # "HBA2","HBB","HBD"
)

missing_genes <- setdiff(select_markers_genes, rownames(samples))
if (length(missing_genes) > 0) {
  warning("Not in object, dropping: ", paste(missing_genes, collapse = ", "))
  select_markers_genes <- setdiff(select_markers_genes, missing_genes)
}
stopifnot(length(select_markers_genes) > 0)

exclude_types <- c("Bipolar TEC","Endothelial","Pericyte","Erythroid")

samples_select_subset <- subset(
  samples,
  subset = !(cell_type_unintegrated %in% exclude_types)
)
samples_select_subset$cell_type_unintegrated <- droplevels(samples_select_subset$cell_type_unintegrated)

# Sanity check
cat("Cell types remaining:\n")
print(table(samples_select_subset$cell_type_unintegrated))
stopifnot(!any(levels(samples_select_subset$cell_type_unintegrated) %in% exclude_types))

cluster_order_select <- c("Nurse", "cTEC", "BMC TEC", "mTEC lo", "mTEC hi",
                          "Mimetic")
stopifnot(setequal(cluster_order_select,
                   levels(samples_select_subset$cell_type_unintegrated)))

avg_mat_select <- AverageExpression(
  samples_select_subset, features = select_markers_genes, assay = "RNA",
  group.by = "cell_type_unintegrated", verbose = TRUE
)$RNA
avg_mat_select <- avg_mat_select[select_markers_genes, cluster_order_select, drop = FALSE]

log_mat_select <- log1p(avg_mat_select)
z_mat_select   <- t(scale(t(log_mat_select)))
z_mat_select[is.nan(z_mat_select)] <- 0

# Sanity check
stopifnot(nrow(z_mat_select) == length(select_markers_genes))
stopifnot(ncol(z_mat_select) == length(cluster_order_select))
stopifnot(!any(is.na(z_mat_select)))

df_heatmap_select <- as.data.frame(z_mat_select) %>%
  tibble::rownames_to_column("gene") %>%
  tidyr::pivot_longer(-gene, names_to = "cluster", values_to = "z") %>%
  mutate(
    gene    = factor(gene, levels = rev(select_markers_genes)),
    cluster = factor(cluster, levels = cluster_order_select)
  )

p_select_heatmap <- ggplot(df_heatmap_select, aes(x = cluster, y = gene, fill = z)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient2(
    low = "#F7E27A", mid = "grey95", high = "#1B3B6F",
    midpoint = 0, name = "Relative\nexpr. (z)"
  ) +
  # Uncomment to display "mTEC I"/"mTEC II" on the axis while keeping
  # mTEC lo/hi as the underlying data values:
  # scale_x_discrete(labels = c("mTEC lo" = "mTEC I", "mTEC hi" = "mTEC II")) +
  labs(title = "Select Markers") +
  theme_minimal(base_size = 10) +
  theme(
    axis.title  = element_blank(),
    plot.title  = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid  = element_blank()
  )
print(p_select_heatmap)
ggsave(paste0(DATE_PREFIX, "-select-markers-heatmap-unintegrated.pdf"),
       p_select_heatmap, width = 10, height = 6)

# -----------------------------------------------------------------------------
# VlnPlot: full "Select Markers" panel, paginated -- 2 cols x 4 rows per page,
# standard 8 x 11.5 in page size
# -----------------------------------------------------------------------------
select_markers_genes <- c(
  "PTPRC", "PDPN",
  "LY75", "PSMB11", "PRSS16",
  "COL4A5", "COL4A6","COL8A1","ITGB4",
  "CD200", "CLU","AIRE", "FEZF2",
  "CD24","ITGA6", "BCAM","CD74",
  "PECAM1", "PDGFRB","RGS5",
  "HBA2","HBB","HBD"
)

missing_genes <- setdiff(select_markers_genes, rownames(samples))
if (length(missing_genes) > 0) {
  warning("Not in object, dropping: ", paste(missing_genes, collapse = ", "))
  select_markers_genes <- setdiff(select_markers_genes, missing_genes)
}
stopifnot(length(select_markers_genes) > 0)

cat("Current cell_type_unintegrated levels:\n")
print(table(samples$cell_type_unintegrated))

genes_per_page <- 8  # 2 cols x 4 rows
n_pages <- ceiling(length(select_markers_genes) / genes_per_page)

# Sanity check
cat(sprintf("\n%d genes across %d page(s), %d genes on the final page\n",
            length(select_markers_genes), n_pages,
            length(select_markers_genes) - genes_per_page * (n_pages - 1)))

pdf_path <- paste0(DATE_PREFIX, "-select-markers-violin-paginated.pdf")
pdf(pdf_path, width = 8.5, height = 11)

for (i in seq_len(n_pages)) {
  start_idx <- (i - 1) * genes_per_page + 1
  end_idx <- min(i * genes_per_page, length(select_markers_genes))
  genes_this_page <- select_markers_genes[start_idx:end_idx]
  
  cat(sprintf("Page %d: %s\n", i, paste(genes_this_page, collapse = ", ")))
  
  p_page <- VlnPlot(
    samples, features = genes_this_page,
    group.by = "cell_type_unintegrated", pt.size = 0.1, ncol = 2
  ) &
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          axis.title.x = element_blank())
  
  print(p_page)
}

dev.off()
cat("\nSaved:", pdf_path, "\n")

# -----------------------------------------------------------------------------
# VlnPlot: full "Select Markers" panel, paginated (2 cols x 4 rows, 8.5x11 in),
# with subsampled point overlay per gene (violin density from ALL cells,
# points shown for only a downsampled subset per group so as not to overwhelm
# the violin plot)
# -----------------------------------------------------------------------------
select_markers_genes <- c(
  "PTPRC", "PDPN",
  "LY75", "PSMB11", "PRSS16",
  "COL4A5", "COL4A6","COL8A1","ITGB4",
  "CD200", "CLU","AIRE", "FEZF2",
  "CD24","ITGA6", "BCAM","CD74",
  "PECAM1", "PDGFRB","RGS5",
  "HBA2","HBB","HBD"
)

group_col <- "cell_type_unintegrated"
n_points_per_group <- 2000   # cap on dots shown per group, per gene

missing_genes <- setdiff(select_markers_genes, rownames(samples))
if (length(missing_genes) > 0) {
  warning("Not in object, dropping: ", paste(missing_genes, collapse = ", "))
  select_markers_genes <- setdiff(select_markers_genes, missing_genes)
}
stopifnot(length(select_markers_genes) > 0)

cat("Current cell_type_unintegrated levels:\n")
print(table(samples[[group_col]]))

Idents(samples) <- group_col
stopifnot(identical(as.character(Idents(samples)), 
                    as.character(samples[[group_col]][[1]])))

# Downsampled cells (same subset used as the point layer for every gene --
# consistent across panels, not re-drawn per gene)
downsampled_cells <- WhichCells(samples, downsample = n_points_per_group)
cat("\nPoints per group after downsampling (applies to every gene panel):\n")
print(table(samples[[group_col]][downsampled_cells, ]))

genes_per_page <- 8  # 2 cols x 4 rows
n_pages <- ceiling(length(select_markers_genes) / genes_per_page)

# Sanity check
cat(sprintf("\n%d genes across %d page(s), %d genes on the final page\n",
            length(select_markers_genes), n_pages,
            length(select_markers_genes) - genes_per_page * (n_pages - 1)))

pdf_path <- paste0(DATE_PREFIX, "-select-markers-violin-paginated-sparse-points.pdf")
pdf(pdf_path, width = 8.5, height = 11)

for (i in seq_len(n_pages)) {
  start_idx <- (i - 1) * genes_per_page + 1
  end_idx <- min(i * genes_per_page, length(select_markers_genes))
  genes_this_page <- select_markers_genes[start_idx:end_idx]
  
  cat(sprintf("Page %d: %s\n", i, paste(genes_this_page, collapse = ", ")))
  
  page_plots <- lapply(genes_this_page, function(g) {
    p_base <- VlnPlot(samples, features = g, group.by = group_col, pt.size = 0) +
      NoLegend() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
            axis.title.x = element_blank())
    
    pt_df <- FetchData(samples, vars = c(g, group_col), cells = downsampled_cells)
    colnames(pt_df) <- c("expr", "group")
    
    p_base +
      geom_jitter(data = pt_df, aes(x = group, y = expr),
                  width = 0.3, size = 0.4, alpha = 0.5, inherit.aes = FALSE)
  })
  
  print(wrap_plots(page_plots, ncol = 2))
}

dev.off()
cat("\nSaved:", normalizePath(pdf_path), "\n")

# -----------------------------------------------------------------------------
# FeaturePlots: full "Select Markers" panel, paginated (3 cols x 4 rows,
# 8.5x11 in per page), scCustomize + viridis_light_high, same rasterization
# approach as another script
# -----------------------------------------------------------------------------
select_markers_genes <- c(
  "PTPRC", "PDPN",
  "LY75", "PSMB11", "PRSS16",
  "COL4A5", "COL4A6","COL8A1","ITGB4",
  "CD200", "CLU","AIRE", "FEZF2",
  "CD24","ITGA6", "BCAM","CD74",
  "PECAM1", "PDGFRB","RGS5",
  "HBA2","HBB","HBD"
)

missing_genes <- setdiff(select_markers_genes, rownames(samples))
if (length(missing_genes) > 0) {
  warning("Not in object, dropping: ", paste(missing_genes, collapse = ", "))
  select_markers_genes <- setdiff(select_markers_genes, missing_genes)
}
stopifnot(length(select_markers_genes) > 0)

# -----------------------------------------------------------------------------
# Panel-building function
# -----------------------------------------------------------------------------
make_panel <- function(gene) {
  p <- FeaturePlot_scCustom(
    samples,
    features = gene,
    colors_use = viridis_light_high,
    order       = TRUE,        # plot high-expression cells on top
    pt.size     = 0.1,
    raster      = FALSE,       # we rasterize ourselves via ggrastr for control
    combine     = TRUE
  ) +
    labs(title = gene) +
    theme_void(base_size = 9) +
    theme(
      plot.title       = element_text(face = "italic", size = 10, hjust = 0.5),
      legend.position  = "right",
      legend.key.width = unit(0.3, "cm"),
      legend.key.height = unit(0.6, "cm"),
      legend.text      = element_text(size = 7),
      legend.title     = element_blank()
    )
  
  # Rasterize just the point geom -- text/axes/legend stay vector
  rasterise(p, layers = "Point", dpi = 300)
}

panels <- lapply(select_markers_genes, make_panel)
names(panels) <- select_markers_genes

# Sanity check
stopifnot(length(panels) == length(select_markers_genes))
cat("Panels built:", length(panels), "\n")

# -----------------------------------------------------------------------------
# Paginate: 3 cols x 4 rows = 12 panels per page, 8.5 x 11 in
# -----------------------------------------------------------------------------
NCOL <- 3
NROW <- 4
panels_per_page <- NCOL * NROW
n_pages <- ceiling(length(panels) / panels_per_page)

# Sanity check
cat(sprintf("%d genes across %d page(s), %d genes on the final page\n",
            length(panels), n_pages,
            length(panels) - panels_per_page * (n_pages - 1)))

pdf_path <- paste0(DATE_PREFIX, "-select-markers-featureplots-paginated.pdf")
pdf(pdf_path, width = 8.5, height = 11)

for (i in seq_len(n_pages)) {
  start_idx <- (i - 1) * panels_per_page + 1
  end_idx <- min(i * panels_per_page, length(panels))
  panels_this_page <- panels[start_idx:end_idx]
  
  cat(sprintf("Page %d: %s\n", i, paste(names(panels_this_page), collapse = ", ")))
  
  composite_page <- wrap_plots(panels_this_page, ncol = NCOL)
  print(composite_page)
  
  # Also save each page as its own PNG (multi-page PNG isn't possible the
  # way multi-page PDF is)
  ggsave(paste0(DATE_PREFIX, "-select-markers-featureplots-page", i, ".png"),
         composite_page, width = 8.5, height = 11, units = "in", dpi = 300)
}

dev.off()
cat("\nSaved:", normalizePath(pdf_path), "\n")
































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
# Cell type display order
# -----------------------------------------------------------------------------
cluster_order_unint <- c("Nurse", "cTEC", "BMC TEC", "mTEC lo", "mTEC hi",
                         "Mimetic","Bipolar TEC","Endothelial",
                         "Pericyte","Erythroid")

stopifnot(setequal(cluster_order_unint, levels(samples$cell_type_unintegrated)))
samples$cell_type_unintegrated <- factor(samples$cell_type_unintegrated,
                                         levels = cluster_order_unint)

# -----------------------------------------------------------------------------
# Ragazzini PolyKRT module score -- VlnPlot with subsampled point overlay,
# NO cell type exclusion (violin density from ALL cells, points downsampled)
# -----------------------------------------------------------------------------
stopifnot("PolyKRT_module_score" %in% colnames(samples@meta.data))

group_col <- "cell_type_unintegrated"
n_points_per_group <- 2000   # same cap used for the Select Markers violin panel

cat("Cell types present (no exclusion):\n")
print(table(samples[[group_col]]))

Idents(samples) <- group_col
stopifnot(identical(as.character(Idents(samples)), as.character(samples[[group_col]][[1]])))

downsampled_cells <- WhichCells(samples, downsample = n_points_per_group)
cat("\nPoints per group after downsampling:\n")
print(table(samples[[group_col]][downsampled_cells, ]))

p_base_polykrt <- VlnPlot(samples, features = "PolyKRT_module_score",
                          group.by = group_col, pt.size = 0) +
  labs(title = "Ragazzini PolyKRT module score by cell type (unintegrated, all types)",
       y = "Module score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pt_df_polykrt <- FetchData(samples, vars = c("PolyKRT_module_score", group_col),
                           cells = downsampled_cells)
colnames(pt_df_polykrt) <- c("expr", "group")

p_vln_all <- p_base_polykrt +
  geom_jitter(data = pt_df_polykrt, aes(x = group, y = expr),
              width = 0.3, size = 0.5, alpha = 0.5, inherit.aes = FALSE)

print(p_vln_all)
ggsave(paste0(DATE_PREFIX, "-PolyKRT-module-score-violin-unintegrated-all-sparse-points.pdf"),
       p_vln_all, width = 9, height = 5)

# -----------------------------------------------------------------------------
# Ragazzini PolyKRT module heatmap -- NO cell type exclusion
# -----------------------------------------------------------------------------
# NOTE: cluster order below is auto-detected from current factor levels, since
# cluster 21's relabel (Fibroblast -> cTEC w/ mesenchymal signature) means a
# hardcoded order from the old script would be stale. Review/reorder manually
# if the printed order below isn't what you want for the figure.
cluster_order_all <- levels(droplevels(samples$cell_type_unintegrated))
cat("Auto-detected cell type order (edit cluster_order_all manually if needed):\n")
print(cluster_order_all)

stopifnot(setequal(cluster_order_all, levels(droplevels(samples$cell_type_unintegrated))))

avg_mat_all <- AverageExpression(
  samples, features = polyKRT_module_genes,
  assay = "RNA", group.by = "cell_type_unintegrated", verbose = TRUE
)$RNA
avg_mat_all <- avg_mat_all[polyKRT_module_genes, cluster_order_all, drop = FALSE]

log_mat_all <- log1p(avg_mat_all)
z_mat_all   <- t(scale(t(log_mat_all)))
z_mat_all[is.nan(z_mat_all)] <- 0

score_by_cluster_all <- samples@meta.data %>%
  group_by(cell_type_unintegrated) %>%
  summarise(mean_score = mean(PolyKRT_module_score), .groups = "drop") %>%
  tibble::deframe()
score_by_cluster_all <- score_by_cluster_all[cluster_order_all]

score_z_all <- as.numeric(scale(score_by_cluster_all))
score_z_all[is.nan(score_z_all)] <- 0
names(score_z_all) <- cluster_order_all

z_mat_full_all <- rbind(z_mat_all, PolyKRT_module_score = score_z_all)
gene_order_full_all <- c(polyKRT_module_genes, "PolyKRT_module_score")

# Sanity check
stopifnot(nrow(z_mat_full_all) == length(gene_order_full_all))
stopifnot(ncol(z_mat_full_all) == length(cluster_order_all))
stopifnot(!any(is.na(z_mat_full_all)))

df_heatmap_all <- as.data.frame(z_mat_full_all) %>%
  tibble::rownames_to_column("gene") %>%
  tidyr::pivot_longer(-gene, names_to = "cluster", values_to = "z") %>%
  mutate(
    gene    = factor(gene, levels = rev(gene_order_full_all)),
    cluster = factor(cluster, levels = cluster_order_all)
  )

p_heatmap_all <- ggplot(df_heatmap_all, aes(x = cluster, y = gene, fill = z)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient2(
    low = "#F7E27A", mid = "grey95", high = "#1B3B6F",
    midpoint = 0, name = "Relative\nexpr. (z)"
  ) +
  labs(title = "PolyKRT module genes + composite score, by cell type (unintegrated, all types)") +
  theme_minimal(base_size = 10) +
  theme(
    axis.title  = element_blank(),
    plot.title  = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid  = element_blank()
  )
print(p_heatmap_all)
ggsave(paste0(DATE_PREFIX, "-PolyKRT-module-heatmap-unintegrated-all.pdf"),
       p_heatmap_all, width = 10, height = 6)

# -----------------------------------------------------------------------------
# Ragazzini PolyKRT module score FeaturePlot
# -----------------------------------------------------------------------------
p_umap <- FeaturePlot(samples, features = "PolyKRT_module_score",
                      reduction = "umap_unintegrated", order = TRUE) +
  scale_color_viridis_c(option = "viridis") +
  labs(title = "PolyKRT module score, UMAP (unintegrated)")
print(p_umap)
ggsave(paste0(DATE_PREFIX, "-PolyKRT-module-score-umap-unintegrated.pdf"),
       p_umap, width = 6, height = 5)

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
cat("  a single population here.\n")





























# -----------------------------------------------------------------------------
# Bautista immature TEC marker heatmap -- unintegrated, subset
# (Bipolar TEC, Pericyte, Fibroblast, Erythroid excluded)
# -----------------------------------------------------------------------------
#exclude_types <- c("Endothelial", "Pericyte", "Erythroid")
exclude_types <- NA

bautista_genes <- c("IGFBP6","KRT15","ZBED2","CDH13","CCN2","DLK2","IGFBP5",
                    "NNMT","MAOA","DPYS","FKBP5","GLUL")

samples_vln_subset <- subset(
  samples,
  subset = !(cell_type_unintegrated %in% exclude_types)
)
samples_vln_subset$cell_type_unintegrated <- droplevels(samples_vln_subset$cell_type_unintegrated)

cat("Cell types remaining:\n")
print(table(samples_vln_subset$cell_type_unintegrated))
stopifnot(!any(levels(samples_vln_subset$cell_type_unintegrated) %in% exclude_types))

# update cell types in cluster order subset below accordingly based on if
# any cell types are excluded in the 'exclude_types' object above
cluster_order_subset <- c("Nurse", "cTEC", "BMC TEC", "mTEC lo", "mTEC hi",
                          "Mimetic", "Bipolar TEC","Endothelial","Pericyte",
                          "Erythroid")

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
# Nurse subcluster characterization -- clusters 3, 13, 15
# Marker enrichment (already computed in markers_unint) + cell cycle scoring
# -----------------------------------------------------------------------------
markers_unint <- read.csv(paste0("results/260807-1222-unintegrated-markers-FULL.csv"))

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
# Nurse subcluster (3, 13, 15) marker analysis: individual PAIRWISE comparisons
# (not pooled) -- 3 vs 13, 3 vs 15, 13 vs 15, each direction separately,
# matching the only.pos=TRUE convention used for the c7 sub-region analysis
# -----------------------------------------------------------------------------
nurse_clusters <- c("3", "13", "15")

samples_nurse_only <- subset(samples, subset = clusters_unintegrated %in% nurse_clusters)
samples_nurse_only$clusters_unintegrated <- droplevels(samples_nurse_only$clusters_unintegrated)

# Sanity check
cat("Cells per Nurse cluster:\n")
print(table(samples_nurse_only$clusters_unintegrated))
stopifnot(setequal(levels(samples_nurse_only$clusters_unintegrated), nurse_clusters))

Idents(samples_nurse_only) <- "clusters_unintegrated"
stopifnot(identical(as.character(Idents(samples_nurse_only)),
                    as.character(samples_nurse_only$clusters_unintegrated)))

# All directed pairs: (3v13, 13v3), (3v15, 15v3), (13v15, 15v13)
pairs <- list(
  c("3", "13"), c("13", "3"),
  c("3", "15"), c("15", "3"),
  c("13", "15"), c("15", "13")
)

results_nurse_pairwise <- list()

for (p in pairs) {
  ident1 <- p[1]; ident2 <- p[2]
  label <- paste0(ident1, "-vs-", ident2)
  
  cat("\n", strrep("=", 60), "\n")
  cat("NURSE CLUSTER", ident1, "vs. CLUSTER", ident2, "(pairwise, not pooled)\n")
  cat(strrep("=", 60), "\n")
  
  n1 <- sum(samples_nurse_only$clusters_unintegrated == ident1)
  n2 <- sum(samples_nurse_only$clusters_unintegrated == ident2)
  cat("Cell counts:", ident1, "=", n1, "|", ident2, "=", n2, "\n")
  
  m <- tryCatch({
    FindMarkers(samples_nurse_only, ident.1 = ident1, ident.2 = ident2,
                only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = TRUE)
  }, error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); NULL })
  
  if (!is.null(m)) {
    m$gene <- rownames(m)
    m <- m %>% arrange(desc(avg_log2FC))
    results_nurse_pairwise[[label]] <- m
    cat("\nTop 15, up in", ident1, "relative to", ident2, ":\n")
    print(head(m[, c("gene","avg_log2FC","pct.1","pct.2","p_val_adj")], 15))
    write.csv(m, paste0(DATE_PREFIX, "-nurse-", label, "-pairwise.csv"),
              row.names = FALSE)
  }
}

# Sanity check
stopifnot(length(results_nurse_pairwise) == length(pairs))
cat("\nAll", length(pairs), "directed pairwise comparisons completed.\n")

# -----------------------------------------------------------------------------
# VlnPlot: Nurse clusters 3, 13, 15 -- conventional layout (1 panel per gene,
# clusters on x-axis within each panel), 3 cols x 3 rows, with subsampled
# point overlay
# -----------------------------------------------------------------------------
genes_to_plot <- c("PRSS16", "PSMB11", "TBATA", "LY75", "RAG1",
                   "CD3E", "MKI67", "TOP2A")
nurse_clusters <- c("3", "13", "15")
n_points_per_group <- 1500

missing_genes <- setdiff(genes_to_plot, rownames(samples))
if (length(missing_genes) > 0) {
  warning("Not in object, dropping: ", paste(missing_genes, collapse = ", "))
  genes_to_plot <- setdiff(genes_to_plot, missing_genes)
}
stopifnot(length(genes_to_plot) > 0)

samples_nurse_only <- subset(samples, subset = clusters_unintegrated %in% nurse_clusters)
samples_nurse_only$clusters_unintegrated <- droplevels(samples_nurse_only$clusters_unintegrated)
samples_nurse_only$clusters_unintegrated <- factor(
  samples_nurse_only$clusters_unintegrated, levels = nurse_clusters
)

# Sanity check
cat("Cells per Nurse cluster:\n")
print(table(samples_nurse_only$clusters_unintegrated))
stopifnot(setequal(levels(samples_nurse_only$clusters_unintegrated), nurse_clusters))

Idents(samples_nurse_only) <- "clusters_unintegrated"
stopifnot(identical(as.character(Idents(samples_nurse_only)),
                    as.character(samples_nurse_only$clusters_unintegrated)))

downsampled_cells <- WhichCells(samples_nurse_only, downsample = n_points_per_group)
cat("\nPoints per cluster after downsampling:\n")
print(table(samples_nurse_only$clusters_unintegrated[downsampled_cells]))

# -----------------------------------------------------------------------------
# Build one panel per gene: violin from all cells, jittered points from
# downsampled subset
# -----------------------------------------------------------------------------
panels <- lapply(genes_to_plot, function(g) {
  p_base <- VlnPlot(samples_nurse_only, features = g,
                    group.by = "clusters_unintegrated", pt.size = 0) +
    NoLegend() +
    labs(x = "Nurse cluster") +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8),
          axis.title.x = element_blank())
  
  pt_df <- FetchData(samples_nurse_only, vars = c(g, "clusters_unintegrated"),
                     cells = downsampled_cells)
  colnames(pt_df) <- c("expr", "group")
  
  p_base +
    geom_jitter(data = pt_df, aes(x = group, y = expr),
                width = 0.3, size = 0.4, alpha = 0.5, inherit.aes = FALSE)
})
names(panels) <- genes_to_plot

# Sanity check
stopifnot(length(panels) == length(genes_to_plot))

p_grid <- wrap_plots(panels, ncol = 4, nrow = 2)
print(p_grid)
ggsave(paste0(DATE_PREFIX, "-nurse-3-13-15-conventional-violin-3x3-sparse-points.pdf"),
       p_grid, width = 8, height = 4)

# -----------------------------------------------------------------------------
# VlnPlot: Nurse clusters 3, 13, 15 -- 1 row per cluster, genes along x-axis,
# with subsampled point overlay (violin density from ALL cells, points
# downsampled per cluster)
# -----------------------------------------------------------------------------
genes_to_plot <- c("RAG1", "CD3E", "MKI67", "TOP2A", "CENPB", "FLT3", "BATF3", 
                   "IGKC", "MS4A1")
nurse_clusters <- c("3", "13", "15")
n_points_per_group <- 1500   # same subsampling cap used before

missing_genes <- setdiff(genes_to_plot, rownames(samples))
if (length(missing_genes) > 0) {
  warning("Not in object, dropping: ", paste(missing_genes, collapse = ", "))
  genes_to_plot <- setdiff(genes_to_plot, missing_genes)
}
stopifnot(length(genes_to_plot) > 0)

samples_nurse_only <- subset(samples, subset = clusters_unintegrated %in% nurse_clusters)
samples_nurse_only$clusters_unintegrated <- droplevels(samples_nurse_only$clusters_unintegrated)
samples_nurse_only$clusters_unintegrated <- factor(
  samples_nurse_only$clusters_unintegrated, levels = nurse_clusters
)

# Sanity check
cat("Cells per Nurse cluster:\n")
print(table(samples_nurse_only$clusters_unintegrated))
stopifnot(setequal(levels(samples_nurse_only$clusters_unintegrated), nurse_clusters))

# -----------------------------------------------------------------------------
# Full data (all cells) -- drives the violin shapes
# -----------------------------------------------------------------------------
expr_data_full <- FetchData(samples_nurse_only,
                            vars = c(genes_to_plot, "clusters_unintegrated"))

df_long_full <- expr_data_full %>%
  tidyr::pivot_longer(cols = all_of(genes_to_plot),
                      names_to = "gene", values_to = "expression") %>%
  mutate(
    gene = factor(gene, levels = genes_to_plot),
    clusters_unintegrated = factor(clusters_unintegrated, levels = nurse_clusters)
  )

# Sanity check
stopifnot(nrow(df_long_full) == ncol(samples_nurse_only) * length(genes_to_plot))
stopifnot(!any(is.na(df_long_full$expression)))

# -----------------------------------------------------------------------------
# Downsampled cells -- drives the point overlay only
# -----------------------------------------------------------------------------
Idents(samples_nurse_only) <- "clusters_unintegrated"
stopifnot(identical(as.character(Idents(samples_nurse_only)),
                    as.character(samples_nurse_only$clusters_unintegrated)))

downsampled_cells <- WhichCells(samples_nurse_only, downsample = n_points_per_group)
cat("\nPoints per cluster after downsampling:\n")
print(table(samples_nurse_only$clusters_unintegrated[downsampled_cells]))

expr_data_sub <- FetchData(samples_nurse_only,
                           vars = c(genes_to_plot, "clusters_unintegrated"),
                           cells = downsampled_cells)

df_long_sub <- expr_data_sub %>%
  tidyr::pivot_longer(cols = all_of(genes_to_plot),
                      names_to = "gene", values_to = "expression") %>%
  mutate(
    gene = factor(gene, levels = genes_to_plot),
    clusters_unintegrated = factor(clusters_unintegrated, levels = nurse_clusters)
  )

# -----------------------------------------------------------------------------
# Plot: violin from full data, jittered points from downsampled subset
# -----------------------------------------------------------------------------
p_nurse_genes_by_row <- ggplot(df_long_full, aes(x = gene, y = expression, fill = clusters_unintegrated)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_jitter(data = df_long_sub, aes(x = gene, y = expression),
              width = 0.15, size = 0.4, alpha = 0.5, color = "black",
              inherit.aes = FALSE) +
  facet_grid(rows = vars(clusters_unintegrated), scales = "free_y") +
  labs(title = "Nurse clusters 3/13/15: select gene expression",
       x = "Gene", y = "Expression level") +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.y = element_text(angle = 0, face = "bold")
  )

print(p_nurse_genes_by_row)
ggsave(paste0(DATE_PREFIX, "-nurse-3-13-15-genes-by-row-violin-sparse-points.pdf"),
       p_nurse_genes_by_row, width = 9, height = 8)

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
# VlnPlot: other markers, Nurse clusters only (3, 13, 15),
# paginated (4 cols x 5 rows, 8.5x11 in), with subsampled point overlay
# -----------------------------------------------------------------------------
select_markers_genes <- c(
  "CCR9", "CD3E", "CD4", "CD8A", "RAG1", "RAG2", "THEMIS",
  "CD1A", "CD1B", "CD1C", "CD1E", "CD38", "BCL11B",
  "DNTT", "PTCRA", "TRBC1", "TRBC2", "SH2D1A",
  "MPEG1", "FLT3", "CCR7", "LAMP3", "BATF3",
  "MS4A1", "CD79A", "CD79B", "PRDM1",
  "LY75", "PSMB11", "PRSS16", "TBATA",
  "EPCAM", "FOXN1"
)

nurse_clusters <- c("3", "13", "15")
group_col <- "clusters_unintegrated"
n_points_per_group <- 1500   # dot-subsampling cap, same as before -- adjust as needed

missing_genes <- setdiff(select_markers_genes, rownames(samples))
if (length(missing_genes) > 0) {
  warning("Not in object, dropping: ", paste(missing_genes, collapse = ", "))
  select_markers_genes <- setdiff(select_markers_genes, missing_genes)
}
stopifnot(length(select_markers_genes) > 0)

samples_nurse_only <- subset(samples, subset = clusters_unintegrated %in% nurse_clusters)
samples_nurse_only$clusters_unintegrated <- droplevels(samples_nurse_only$clusters_unintegrated)
samples_nurse_only$clusters_unintegrated <- factor(
  samples_nurse_only$clusters_unintegrated, levels = nurse_clusters
)

# Sanity check
cat("Cells per Nurse cluster:\n")
print(table(samples_nurse_only$clusters_unintegrated))
stopifnot(setequal(levels(samples_nurse_only$clusters_unintegrated), nurse_clusters))

Idents(samples_nurse_only) <- group_col
stopifnot(identical(as.character(Idents(samples_nurse_only)),
                    as.character(samples_nurse_only[[group_col]][[1]])))

downsampled_cells <- WhichCells(samples_nurse_only, downsample = n_points_per_group)
cat("\nPoints per group after downsampling:\n")
print(table(samples_nurse_only[[group_col]][downsampled_cells, ]))

genes_per_page <- 20  # 4 cols x 5 rows
n_pages <- ceiling(length(select_markers_genes) / genes_per_page)

# Sanity check
cat(sprintf("\n%d genes across %d page(s), %d genes on the final page\n",
            length(select_markers_genes), n_pages,
            length(select_markers_genes) - genes_per_page * (n_pages - 1)))

pdf_path <- paste0(DATE_PREFIX, "-nurse-clusters-select-markers-violin-paginated.pdf")
pdf(pdf_path, width = 8.5, height = 11)

for (i in seq_len(n_pages)) {
  start_idx <- (i - 1) * genes_per_page + 1
  end_idx <- min(i * genes_per_page, length(select_markers_genes))
  genes_this_page <- select_markers_genes[start_idx:end_idx]
  
  cat(sprintf("Page %d: %s\n", i, paste(genes_this_page, collapse = ", ")))
  
  page_plots <- lapply(genes_this_page, function(g) {
    p_base <- VlnPlot(samples_nurse_only, features = g,
                      group.by = group_col, pt.size = 0) +
      NoLegend() +
      labs(x = "Nurse cluster") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
            axis.title.x = element_blank())
    
    pt_df <- FetchData(samples_nurse_only, vars = c(g, group_col),
                       cells = downsampled_cells)
    colnames(pt_df) <- c("expr", "group")
    
    p_base +
      geom_jitter(data = pt_df, aes(x = group, y = expr),
                  width = 0.3, size = 0.4, alpha = 0.5, inherit.aes = FALSE)
  })
  
  print(wrap_plots(page_plots, ncol = 4))
}

dev.off()
cat("\nSaved:", normalizePath(pdf_path), "\n")

# =============================================================================
# Cluster 15 standalone re-embedding -- same method as Bipolar TEC and Nurse
# subclustering: subset, wipe stale reductions, fresh PCA/UMAP/clustering on
# just these cells.
# =============================================================================

# -----------------------------------------------------------------------------
# Step 1: Subset to cluster 15 cells only
# -----------------------------------------------------------------------------
cluster15_standalone <- subset(samples, subset = clusters_unintegrated == "15")

cat("Cluster 15 standalone cell count:", ncol(cluster15_standalone), "\n")

# -----------------------------------------------------------------------------
# Step 2: Wipe stale layers/reductions from the parent object
# -----------------------------------------------------------------------------
DefaultAssay(cluster15_standalone) <- "RNA"
cluster15_standalone <- DietSeurat(
  cluster15_standalone,
  layers    = c("counts", "data"),
  dimreducs = NULL,
  graphs    = NULL,
  misc      = FALSE
)

# Sanity check
stopifnot(length(Reductions(cluster15_standalone)) == 0)
cat("\nLayers after wipe:", paste(Layers(cluster15_standalone, assay = "RNA"), collapse = ", "), "\n")

# -----------------------------------------------------------------------------
# Step 3: Fresh preprocessing on just these cells
# -----------------------------------------------------------------------------
cluster15_standalone <- NormalizeData(cluster15_standalone, verbose = TRUE)
cluster15_standalone <- FindVariableFeatures(cluster15_standalone, nfeatures = 2000, verbose = TRUE)
cluster15_standalone <- ScaleData(cluster15_standalone, verbose = TRUE)
cluster15_standalone <- RunPCA(cluster15_standalone, npcs = 30, verbose = TRUE)

p_elbow_c15 <- ElbowPlot(cluster15_standalone, ndims = 30) +
  labs(title = "Cluster 15 standalone PCA: elbow plot")
print(p_elbow_c15)
ggsave(paste0(DATE_PREFIX, "-cluster15-standalone-elbow.pdf"), p_elbow_c15, width = 6, height = 5)

# -----------------------------------------------------------------------------
# Step 4: Neighbors, clustering (sweep a few resolutions), UMAP
# -----------------------------------------------------------------------------
N_PCS <- 10

cluster15_standalone <- FindNeighbors(cluster15_standalone, reduction = "pca",
                                      dims = 1:N_PCS, verbose = TRUE)

resolutions <- c(0.1, 0.2, 0.4, 0.6, 0.8)
for (r in resolutions) {
  cluster15_standalone <- FindClusters(cluster15_standalone, resolution = r, verbose = TRUE)
}

cluster15_standalone <- RunUMAP(cluster15_standalone, reduction = "pca",
                                dims = 1:N_PCS, verbose = TRUE)

cat("\n=== Cluster counts at each resolution ===\n")
for (r in resolutions) {
  col <- paste0("RNA_snn_res.", r)
  cat(sprintf("res %.1f -> %d clusters\n", r, length(unique(cluster15_standalone[[col]][[1]]))))
}

cat("\n=== Cluster sizes at each resolution ===\n")
for (r in resolutions) {
  col <- paste0("RNA_snn_res.", r)
  cat("\n--- res =", r, "---\n")
  print(table(cluster15_standalone[[col]][[1]]))
}

# -----------------------------------------------------------------------------
# Step 5: Cross-tab data-driven clusters (res=0.1 and res=0.2) against
# FLT3+/MS4A1+ status
# -----------------------------------------------------------------------------
fm_status <- FetchData(cluster15_standalone, vars = c("FLT3", "MS4A1")) %>%
  mutate(
    status = case_when(
      FLT3 > 0 & MS4A1 > 0 ~ "double_pos",
      FLT3 > 0 & MS4A1 == 0 ~ "FLT3_pos",
      FLT3 == 0 & MS4A1 > 0 ~ "MS4A1_pos",
      TRUE ~ "double_neg"
    )
  )

# Sanity check
stopifnot(nrow(fm_status) == ncol(cluster15_standalone))
cat("Overall FLT3/MS4A1 status counts:\n")
print(table(fm_status$status))

cluster15_standalone$fm_status <- fm_status$status

for (r in c(0.1, 0.2)) {
  col <- paste0("RNA_snn_res.", r)
  cat("\n=== res =", r, "-- data-driven cluster vs. FLT3/MS4A1 status ===\n")
  cross_tab <- table(cluster15_standalone[[col]][[1]], cluster15_standalone$fm_status)
  print(cross_tab)
  
  # Sanity check
  stopifnot(sum(cross_tab) == ncol(cluster15_standalone))
  
  # Row percentages -- purity of each data-driven cluster
  cat("\nRow percentages (purity within each data-driven cluster):\n")
  print(round(100 * prop.table(cross_tab, margin = 1), 1))
}

# -----------------------------------------------------------------------------
# Replace single-gene FLT3/MS4A1 classification with multi-gene module scores
# -----------------------------------------------------------------------------
dc_genes    <- c("FLT3", "BATF3", "WDFY4", "CCR7", "LAMP3", "FSCN1")
tolero_genes <- c("CD80", "CD86", "IDO1", "KYNU", "IL4I1")
pdc_genes   <- c("SPIB", "CLEC4C", "IL3RA", "LILRA4", "TCF4", "IRF7")
b_genes     <- c("MS4A1", "CD79A", "CD79B", "IGKC", "BTK")
myeloid_genes <- c("MPEG1", "LYZ", "CYBB", "CPVL", "FGL2")

module_list <- list(dc_genes, tolero_genes, pdc_genes, b_genes, myeloid_genes)
module_names <- c("DC", "Tolerogenic", "pDC", "Bcell", "Myeloid")

# sanity check: confirm all module genes are present in the object before scoring
for (i in seq_along(module_list)) {
  missing <- setdiff(module_list[[i]], rownames(cluster15_standalone))
  if (length(missing) > 0) {
    warning(module_names[i], " module missing genes: ", paste(missing, collapse = ", "))
    module_list[[i]] <- setdiff(module_list[[i]], missing)
  }
}

cluster15_standalone <- AddModuleScore(
  cluster15_standalone, features = module_list, name = "LineageScore", verbose = TRUE
)
# AddModuleScore appends 1:length(module_list) as suffix -- rename for clarity
for (i in seq_along(module_names)) {
  colnames(cluster15_standalone@meta.data)[
    colnames(cluster15_standalone@meta.data) == paste0("LineageScore", i)
  ] <- paste0(module_names[i], "_score")
}

# sanity check: confirm all 5 module score columns were created with no NAs
score_cols <- paste0(module_names, "_score")
stopifnot(all(score_cols %in% colnames(cluster15_standalone@meta.data)))
stopifnot(!any(is.na(cluster15_standalone@meta.data[, score_cols])))

# Assign each cell to its highest-scoring lineage
lineage_calls <- apply(cluster15_standalone@meta.data[, score_cols], 1, function(x) {
  module_names[which.max(x)]
})
cluster15_standalone$lineage_call <- lineage_calls

cat("Lineage assignment counts (by max module score):\n")
print(table(cluster15_standalone$lineage_call))

# Cross-tab against res=0.2 data-driven clusters
cat("\n=== res=0.2 cluster vs. module-score lineage call ===\n")
cross_tab_module <- table(cluster15_standalone$RNA_snn_res.0.2, cluster15_standalone$lineage_call)
print(cross_tab_module)
print(round(100 * prop.table(cross_tab_module, margin = 1), 1))

# -----------------------------------------------------------------------------
# Add cTEC module score -- kept SEPARATE from the immune lineage "max score"
# classification, since cTEC signal (real or contaminating) isn't mutually
# exclusive with any of DC/B/pDC/Myeloid/Tolerogenic identity
# -----------------------------------------------------------------------------
ctec_genes <- c("LY75", "PRSS16", "PSMB11", "TBATA", "CCL25", "DLL4")

missing_ctec <- setdiff(ctec_genes, rownames(cluster15_standalone))
if (length(missing_ctec) > 0) {
  warning("cTEC module missing genes: ", paste(missing_ctec, collapse = ", "))
  ctec_genes <- setdiff(ctec_genes, missing_ctec)
}
stopifnot(length(ctec_genes) > 0)

cluster15_standalone <- AddModuleScore(
  cluster15_standalone, features = list(ctec_genes), name = "cTEC_score_raw", verbose = TRUE
)
cluster15_standalone$cTEC_score <- cluster15_standalone$cTEC_score_raw1
cluster15_standalone$cTEC_score_raw1 <- NULL

# sanity check: confirm cTEC_score column created with no NAs
stopifnot("cTEC_score" %in% colnames(cluster15_standalone@meta.data))
stopifnot(!any(is.na(cluster15_standalone$cTEC_score)))

# -----------------------------------------------------------------------------
# cTEC score by res=0.2 data-driven cluster
# -----------------------------------------------------------------------------
cat("=== cTEC module score by res=0.2 cluster ===\n")
ctec_by_res02 <- cluster15_standalone@meta.data %>%
  group_by(RNA_snn_res.0.2) %>%
  summarise(
    n = n(),
    mean_cTEC_score = round(mean(cTEC_score), 3),
    median_cTEC_score = round(median(cTEC_score), 3),
    pct_cTEC_score_above_0 = round(100 * mean(cTEC_score > 0), 1),
    .groups = "drop"
  )
print(ctec_by_res02)

# -----------------------------------------------------------------------------
# cTEC score by immune lineage call
# -----------------------------------------------------------------------------
cat("\n=== cTEC module score by immune lineage call (max-score classification) ===\n")
ctec_by_lineage <- cluster15_standalone@meta.data %>%
  group_by(lineage_call) %>%
  summarise(
    n = n(),
    mean_cTEC_score = round(mean(cTEC_score), 3),
    median_cTEC_score = round(median(cTEC_score), 3),
    pct_cTEC_score_above_0 = round(100 * mean(cTEC_score > 0), 1),
    .groups = "drop"
  )
print(ctec_by_lineage)

# -----------------------------------------------------------------------------
# Visualize: cTEC score distribution across res=0.2 clusters
# -----------------------------------------------------------------------------
p_ctec_by_cluster <- VlnPlot(cluster15_standalone, features = "cTEC_score",
                             group.by = "RNA_snn_res.0.2", pt.size = 0.2) +
  labs(title = "cTEC module score across cluster 15 standalone sub-clusters (res=0.2)",
       x = "Sub-cluster", y = "cTEC module score")
print(p_ctec_by_cluster)
ggsave(paste0(DATE_PREFIX, "-cluster15-cTEC-score-by-subcluster.pdf"),
       p_ctec_by_cluster, width = 4, height = 3)

# Pairwise Wilcoxon: does any subcluster's cTEC score differ significantly
# from the rest?
cat("\n=== Wilcoxon test: each res=0.2 cluster's cTEC score vs. all others ===\n")
for (cl in levels(cluster15_standalone$RNA_snn_res.0.2)) {
  in_cluster <- cluster15_standalone$cTEC_score[cluster15_standalone$RNA_snn_res.0.2 == cl]
  out_cluster <- cluster15_standalone$cTEC_score[cluster15_standalone$RNA_snn_res.0.2 != cl]
  test <- wilcox.test(in_cluster, out_cluster)
  cat(sprintf("Cluster %s: mean=%.3f vs rest mean=%.3f, p=%.2e\n",
              cl, mean(in_cluster), mean(out_cluster), test$p.value))
}

# -----------------------------------------------------------------------------
# VlnPlot: individual cTEC genes across cluster 15 sub-clusters (res=0.2)
# -----------------------------------------------------------------------------
ctec_genes_check <- c("LY75", "PRSS16", "PSMB11", "TBATA", "CCL25", "DLL4")

missing_ctec <- setdiff(ctec_genes_check, rownames(cluster15_standalone))
if (length(missing_ctec) > 0) {
  warning("Not in object, dropping: ", paste(missing_ctec, collapse = ", "))
  ctec_genes_check <- setdiff(ctec_genes_check, missing_ctec)
}
stopifnot(length(ctec_genes_check) > 0)

p_ctec_genes_c15 <- VlnPlot(
  cluster15_standalone, features = ctec_genes_check,
  group.by = "RNA_snn_res.0.2", pt.size = 0.2, ncol = 3
) &
  labs(x = "Sub-cluster")

print(p_ctec_genes_c15)
ggsave(paste0(DATE_PREFIX, "-cluster15-ctec-genes-by-subcluster-violin.pdf"),
       p_ctec_genes_c15, width = 12, height = 8)

# -----------------------------------------------------------------------------
# Vlnplot of essential defining markers for the cluster 15 subclusters
# -----------------------------------------------------------------------------
condensed_genes <- c(
  "PRSS16","TBATA","FLT3","TRAC","RAG1","MKI67",
  "CLEC4C","CD19","IGKC","XCR1","C1QA","CD163"
  # "MS4A1","CD79A","CD79B","IGKC",
  # "FLT3","IRF4","IRF7","IRF8",
  # "SPIB","XCR1","CLEC9A",
  # "GIMAP7","CD3G","THEMIS","CD7"
  # "C1QA","CD163","CSF1R","MKI67"
)
missing_genes <- setdiff(condensed_genes, rownames(cluster15_standalone))
if (length(missing_genes) > 0) {
  warning("Not in object, dropping: ", paste(missing_genes, collapse = ", "))
  condensed_genes <- setdiff(condensed_genes, missing_genes)
}
stopifnot(length(condensed_genes) > 0)
p_condensed <- VlnPlot(
  cluster15_standalone, features = condensed_genes,
  group.by = "RNA_snn_res.0.2", pt.size = 0.2, ncol = 4
) &
  labs(x = "Sub-cluster") &
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
print(p_condensed)
ggsave(paste0(DATE_PREFIX, "-cluster15-condensed-8panel-violin.pdf"),
       p_condensed, width = 8, height = 6)

# -----------------------------------------------------------------------------
# Sub-cluster 2 defining genes -- vs. rest of cluster15_standalone
# -----------------------------------------------------------------------------
Idents(cluster15_standalone) <- "RNA_snn_res.0.2"
# sanity check: confirm Idents reflects res=0.2 clustering before FindMarkers
stopifnot(identical(as.character(Idents(cluster15_standalone)),
                    as.character(cluster15_standalone$RNA_snn_res.0.2)))

markers_sub2 <- FindMarkers(cluster15_standalone, ident.1 = "2", only.pos = TRUE,
                            min.pct = 0.25, logfc.threshold = 0.25, verbose = TRUE)
markers_sub2$gene <- rownames(markers_sub2)
markers_sub2 <- markers_sub2 %>% arrange(desc(avg_log2FC))

cat("Top 20, sub-cluster 2:\n")
print(head(markers_sub2[, c("gene","avg_log2FC","pct.1","pct.2","p_val_adj")], 20))

write.csv(markers_sub2, paste0(DATE_PREFIX, "-cluster15-subcluster2-markers.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# Direct pairwise comparison: cluster 3 vs cluster 4 (cluster15_standalone,
# res=0.2) -- what distinguishes the "clean B cell" cluster from the
# "B cell/pDC hybrid" cluster
# -----------------------------------------------------------------------------
Idents(cluster15_standalone) <- "RNA_snn_res.0.2"
# sanity check: confirm Idents reflects res=0.2 clustering before FindMarkers
stopifnot(identical(as.character(Idents(cluster15_standalone)),
                    as.character(cluster15_standalone$RNA_snn_res.0.2)))

markers_3v4 <- FindMarkers(cluster15_standalone, ident.1 = "3", ident.2 = "4",
                           only.pos = TRUE, min.pct = 0.25,
                           logfc.threshold = 0.25, verbose = TRUE)
markers_3v4$gene <- rownames(markers_3v4)
markers_3v4 <- markers_3v4 %>% arrange(desc(avg_log2FC))

markers_4v3 <- FindMarkers(cluster15_standalone, ident.1 = "4", ident.2 = "3",
                           only.pos = TRUE, min.pct = 0.25,
                           logfc.threshold = 0.25, verbose = TRUE)
markers_4v3$gene <- rownames(markers_4v3)
markers_4v3 <- markers_4v3 %>% arrange(desc(avg_log2FC))

cat("Top 20, up in cluster 3 relative to cluster 4:\n")
print(head(markers_3v4[, c("gene","avg_log2FC","pct.1","pct.2","p_val_adj")], 20))

cat("\nTop 20, up in cluster 4 relative to cluster 3:\n")
print(head(markers_4v3[, c("gene","avg_log2FC","pct.1","pct.2","p_val_adj")], 20))

write.csv(markers_3v4, paste0(DATE_PREFIX, "-cluster15-sub3-vs-sub4.csv"), row.names = FALSE)
write.csv(markers_4v3, paste0(DATE_PREFIX, "-cluster15-sub4-vs-sub3.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# Direct pairwise comparison: cluster 5 vs cluster 6 (cluster15_standalone,
# res=0.2) -- what actually distinguishes these two myeloid sub-clusters
# -----------------------------------------------------------------------------
Idents(cluster15_standalone) <- "RNA_snn_res.0.2"
# sanity check: confirm Idents reflects res=0.2 clustering before FindMarkers
stopifnot(identical(as.character(Idents(cluster15_standalone)),
                    as.character(cluster15_standalone$RNA_snn_res.0.2)))

markers_5v6 <- FindMarkers(cluster15_standalone, ident.1 = "5", ident.2 = "6",
                           only.pos = TRUE, min.pct = 0.25,
                           logfc.threshold = 0.25, verbose = TRUE)
markers_5v6$gene <- rownames(markers_5v6)
markers_5v6 <- markers_5v6 %>% arrange(desc(avg_log2FC))

markers_6v5 <- FindMarkers(cluster15_standalone, ident.1 = "6", ident.2 = "5",
                           only.pos = TRUE, min.pct = 0.25,
                           logfc.threshold = 0.25, verbose = TRUE)
markers_6v5$gene <- rownames(markers_6v5)
markers_6v5 <- markers_6v5 %>% arrange(desc(avg_log2FC))

cat("Top 15, up in cluster 5 relative to cluster 6:\n")
print(head(markers_5v6[, c("gene","avg_log2FC","pct.1","pct.2","p_val_adj")], 15))

cat("\nTop 15, up in cluster 6 relative to cluster 5:\n")
print(head(markers_6v5[, c("gene","avg_log2FC","pct.1","pct.2","p_val_adj")], 15))

write.csv(markers_5v6, paste0(DATE_PREFIX, "-cluster15-sub5-vs-sub6.csv"), row.names = FALSE)
write.csv(markers_6v5, paste0(DATE_PREFIX, "-cluster15-sub6-vs-sub5.csv"), row.names = FALSE)








































# =============================================================================
# Figure: Cluster 7 analysis
# =============================================================================

# -----------------------------------------------------------------------------
# Panel A: Cluster 7 (red) highlighted on the whole-dataset UMAP (unintegrated)
# -----------------------------------------------------------------------------
p_panelA <- DimPlot(
  samples, reduction = "umap_unintegrated",
  cells.highlight = WhichCells(samples, expression = clusters_unintegrated == "7"),
  sizes.highlight = 0.6, pt.size = 0.3,
  split.by = "sample.id"
) +
  scale_color_manual(values = c("grey85", "red"), labels = c("Other", "Bipolar TEC")) +
  labs(title = "Bipolar TEC (cluster 7) on whole-dataset UMAP") +
  theme_void(base_size = 10) +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
        legend.position = "right")

print(p_panelA)
ggsave(paste0(DATE_PREFIX, "-Fig-BipolarTEC-panelA-umap-highlight-split-by-sample.pdf"),
       p_panelA, width = 8, height = 4)

# sanity check: confirm highlighted cell count matches cluster 7's known size
stopifnot(length(WhichCells(samples, expression = clusters_unintegrated == "7")) == 1264)

# -----------------------------------------------------------------------------
# Cluster 7 QC Test: nFeature, nCount, percent.mt
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# QC Test #2 (re-run on current object): nFeature, nCount, percent.mt --
# side-by-side violin, cluster 7 vs. All, with subsampled point
# overlay (violin density from ALL cells, points downsampled per group)
# -----------------------------------------------------------------------------
samples$c7_vs_all <- ifelse(samples$clusters_unintegrated == "7",
                              "Cluster 7", "All")
samples$c7_vs_all <- factor(samples$c7_vs_all,
                              levels = c("Cluster 7", "All"))

# sanity check: confirm group sizes match cluster 7's known size and total dataset
stopifnot(sum(samples$c7_vs_all == "Cluster 7") == 1264)
stopifnot(sum(table(samples$c7_vs_all)) == ncol(samples))

n_points_per_group <- 5000   # same cap used throughout this session

Idents(samples) <- "c7_vs_all"
# sanity check: confirm Idents matches c7_vs_all before downsampling
stopifnot(identical(as.character(Idents(samples)), as.character(samples$c7_vs_all)))

downsampled_cells <- WhichCells(samples, downsample = n_points_per_group)
cat("Points per group after downsampling:\n")
print(table(samples$c7_vs_all[downsampled_cells]))

qc_features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")

panels_qc <- lapply(qc_features, function(feat) {
  p_base <- VlnPlot(samples, features = feat, group.by = "c7_vs_all", pt.size = 0) +
    NoLegend() +
    labs(x = NULL) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
  
  pt_df <- FetchData(samples, vars = c(feat, "c7_vs_all"), cells = downsampled_cells)
  colnames(pt_df) <- c("value", "group")
  
  p_base +
    geom_jitter(data = pt_df, aes(x = group, y = value),
                width = 0.3, size = 0.4, alpha = 0.5, inherit.aes = FALSE)
})

p_qc_c7_vs_all <- wrap_plots(panels_qc, ncol = 3)
print(p_qc_c7_vs_all)
ggsave(paste0(DATE_PREFIX, "-c7-vs-whole-QC-metrics-violin-sparse-points.pdf"),
       p_qc_c7_vs_all, width = 8, height = 6)

# -----------------------------------------------------------------------------
# Cluster 7 QC Test: nFeature vs nCount scatter highlighted on top of 
# the all-dataset distribution
# -----------------------------------------------------------------------------
cells_c7 <- colnames(samples)[samples$clusters_unintegrated == "7"]

# sanity check: confirm highlighted cell count matches cluster 7's known size
stopifnot(length(cells_c7) == 1264)

qc_df_c7 <- samples@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  select(cell, nFeature_RNA, nCount_RNA) %>%
  mutate(
    highlight = ifelse(cell %in% cells_c7, "Cluster 7", "all dataset"),
    highlight = factor(highlight, levels = c("all dataset", "Cluster 7"))
  )

# sanity check: confirm highlighted count matches cells_c7, rest is background
stopifnot(sum(qc_df_c7$highlight == "Cluster 7") == length(cells_c7))

p_qc_c7_scatter <- ggplot() +
  geom_point(data = filter(qc_df_c7, highlight == "all dataset"),
             aes(x = nCount_RNA, y = nFeature_RNA),
             color = "grey80", size = 0.3, alpha = 0.4) +
  geom_point(data = filter(qc_df_c7, highlight == "Cluster 7"),
             aes(x = nCount_RNA, y = nFeature_RNA),
             color = "red", size = 0.6, alpha = 0.6) +
  labs(title = "nFeature vs nCount: all dataset (grey) with cluster 7 highlighted (red)",
       x = "nCount_RNA", y = "nFeature_RNA") +
  theme_minimal(base_size = 11)

print(p_qc_c7_scatter)
ggsave(paste0(DATE_PREFIX, "-c7-nFeature-vs-nCount-scatter.pdf"),
       p_qc_c7_scatter, width = 8, height = 7)

# -----------------------------------------------------------------------------
# Panel B: Heatmap -- cluster 7's top 15 FindMarkers genes (vs. all other
# clusters), shown across all annotated cell types including Bipolar TEC
# -----------------------------------------------------------------------------
bipolar_top15_genes <- top15_unint %>%
  filter(cluster == "7") %>%
  arrange(desc(avg_log2FC)) %>%
  pull(gene)

# sanity check: confirm exactly 15 genes retrieved for cluster 7
stopifnot(length(bipolar_top15_genes) == 15)
cat("Bipolar TEC top 15 genes:\n")
cat(paste(bipolar_top15_genes, collapse = ", "), "\n")

cluster_order_panelB <- c("Nurse", "cTEC", "BMC TEC", "mTEC lo", "mTEC hi",
                          "Mimetic", "Endothelial", "Bipolar TEC", "Pericyte",
                          "Erythroid")

# sanity check: confirm cluster_order_panelB matches current object's cell type levels
current_levels <- levels(droplevels(samples$cell_type_unintegrated))
cat("\nCurrent cell_type_unintegrated levels:\n")
print(current_levels)
stopifnot(setequal(cluster_order_panelB, current_levels))

avg_mat_panelB <- AverageExpression(
  samples, features = bipolar_top15_genes, assay = "RNA",
  group.by = "cell_type_unintegrated", verbose = TRUE
)$RNA
avg_mat_panelB <- avg_mat_panelB[bipolar_top15_genes, cluster_order_panelB, drop = FALSE]

log_mat_panelB <- log1p(avg_mat_panelB)
z_mat_panelB   <- t(scale(t(log_mat_panelB)))
z_mat_panelB[is.nan(z_mat_panelB)] <- 0

# sanity check: confirm heatmap matrix dimensions match gene count x cell type count
stopifnot(nrow(z_mat_panelB) == 15)
stopifnot(ncol(z_mat_panelB) == length(cluster_order_panelB))
stopifnot(!any(is.na(z_mat_panelB)))

df_panelB <- as.data.frame(z_mat_panelB) %>%
  tibble::rownames_to_column("gene") %>%
  tidyr::pivot_longer(-gene, names_to = "cluster", values_to = "z") %>%
  mutate(
    gene    = factor(gene, levels = rev(bipolar_top15_genes)),
    cluster = factor(cluster, levels = cluster_order_panelB)
  )

p_panelB <- ggplot(df_panelB, aes(x = cluster, y = gene, fill = z)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient2(
    low = "#F7E27A", mid = "grey95", high = "#1B3B6F",
    midpoint = 0, name = "Relative\nexpr. (z)"
  ) +
  labs(title = "Bipolar TEC top 15 markers (vs. all other clusters)") +
  theme_minimal(base_size = 10) +
  theme(
    axis.title  = element_blank(),
    plot.title  = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid  = element_blank()
  )

print(p_panelB)
ggsave(paste0(DATE_PREFIX, "-Fig-BipolarTEC-panelB-top15-heatmap.pdf"),
       p_panelB, width = 7, height = 6)

# -----------------------------------------------------------------------------
# Heatmap: top 15 markers defining Bipolar TEC
# -----------------------------------------------------------------------------
bipolar_genes <- top15_unint %>%
  filter(cluster == "7") %>%
  arrange(desc(avg_log2FC)) %>%
  pull(gene)

stopifnot(length(bipolar_genes) == 15)
cat("Bipolar TEC top 15 genes:\n")
cat(paste(bipolar_genes, collapse = ", "), "\n")

#exclude_types <- c("Endothelial", "Pericyte", "Erythroid")   # Bipolar TEC retained
exclude_types <- NA

samples_bipolar_subset <- subset(
  samples,
  subset = !(cell_type_unintegrated %in% exclude_types)
)
samples_bipolar_subset$cell_type_unintegrated <- droplevels(samples_bipolar_subset$cell_type_unintegrated)

cat("\nCell types remaining:\n")
print(table(samples_bipolar_subset$cell_type_unintegrated))
stopifnot(!any(levels(samples_bipolar_subset$cell_type_unintegrated) %in% exclude_types))

cluster_order_bipolar <- c("Nurse", "cTEC", "BMC TEC", "mTEC lo", "mTEC hi",
                           "Mimetic","Endothelial", "Bipolar TEC", "Pericyte",
                           "Erythroid")

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
stopifnot(identical(as.character(Idents(samples)), as.character(samples$ident_for_c7_analysis)))

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

# -----------------------------------------------------------------------------
# Bipolar TEC sub-region analysis: individual PAIRWISE comparisons (not pooled)
# -----------------------------------------------------------------------------
Idents(samples) <- "ident_for_c7_analysis"
stopifnot(all(as.character(Idents(samples)) == unname(samples$ident_for_c7_analysis)))

subregions <- c("near_cTEC", "near_mTEC_lo", "near_mTEC_hi", "near_Nurse")
region_pairs <- combn(subregions, 2, simplify = FALSE)  # 6 unique pairs

results_pairwise <- list()

for (pair in region_pairs) {
  for (direction in list(c(pair[1], pair[2]), c(pair[2], pair[1]))) {
    ident1 <- direction[1]; ident2 <- direction[2]
    label <- paste0(ident1, "-vs-", ident2)
    
    cat("\n", strrep("=", 60), "\n")
    cat(ident1, "vs.", ident2, "(pairwise, not pooled)\n")
    cat(strrep("=", 60), "\n")
    
    n1 <- sum(samples$ident_for_c7_analysis == ident1, na.rm = TRUE)
    n2 <- sum(samples$ident_for_c7_analysis == ident2, na.rm = TRUE)
    cat("Cell counts:", ident1, "=", n1, "|", ident2, "=", n2, "\n")
    
    m <- tryCatch({
      FindMarkers(samples, ident.1 = ident1, ident.2 = ident2,
                  only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = TRUE)
    }, error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); NULL })
    
    if (!is.null(m)) {
      m$gene <- rownames(m)
      m <- m %>% arrange(desc(avg_log2FC))
      results_pairwise[[label]] <- m
      cat("\nTop 15, up in", ident1, "relative to", ident2, ":\n")
      print(head(m[, c("gene","avg_log2FC","pct.1","pct.2","p_val_adj")], 15))
      write.csv(m, paste0(DATE_PREFIX, "-c7-", label, "-pairwise.csv"), row.names = FALSE)
    }
  }
}

# sanity check: confirm all 12 directed pairwise comparisons completed
stopifnot(length(results_pairwise) == length(region_pairs) * 2)
cat("\nAll", length(results_pairwise), "directed pairwise comparisons completed.\n")

# -----------------------------------------------------------------------------
# VlnPlot: top genes from the c7 sub-region pairwise comparisons, across ALL
# annotated cell types in the All -- check whether each gene's
# expression actually tracks with the corresponding neighbor cluster
# -----------------------------------------------------------------------------
genes_to_check <- c(
  # near_cTEC vs near_mTEC_lo (weak internal separation -- check vs real neighbors)
  "TNFRSF17", "NPTX1", "KCNE5", "TP53AIP1", "CA8",
  # near_mTEC_lo vs near_cTEC (reverse direction)
  "CHI3L1", "SLC34A2", "ITGB8", "CXCL10", "MYL9",
  # near_mTEC_hi's consistently top gene, for contrast
  "FEZF2", "DLX5",
  # near_Nurse's consistently top genes, for contrast
  "ARHGAP15", "RHOH", "CD4"
)

missing_genes <- setdiff(genes_to_check, rownames(samples))
if (length(missing_genes) > 0) {
  warning("Not in object, dropping: ", paste(missing_genes, collapse = ", "))
  genes_to_check <- setdiff(genes_to_check, missing_genes)
}
stopifnot(length(genes_to_check) > 0)

# sanity check: confirm cell_type_unintegrated levels before plotting
cat("Current cell_type_unintegrated levels:\n")
print(table(samples$cell_type_unintegrated))

p_genes_whole_dataset <- VlnPlot(
  samples, features = genes_to_check,
  group.by = "cell_type_unintegrated", pt.size = 0, ncol = 3
) &
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.title.x = element_blank())

print(p_genes_whole_dataset)
ggsave(paste0(DATE_PREFIX, "-c7-subregion-genes-whole-dataset-violin.png"),
       p_genes_whole_dataset, width = 12, height = 14)

# -----------------------------------------------------------------------------
# VlnPlot: top 5 genes from the POOLED vs-OTHER3-internal analysis for each
# c7 sub-region, plotted across ALL cell types in the All
# -----------------------------------------------------------------------------
genes_pooled_check <- c(
  "CD4", "SPN", "CD96", "RHOH", "CD247",              # near_Nurse (pooled)
  "TNFRSF17", "HAO2", "SNX22",                        # near_cTEC (pooled)
  "FEZF2", "DLX5", "CLDN3", "SRGN", "CXCL5",         # near_mTEC_hi (pooled)
  "MYL9", "CXCL10", "C3", "SLC34A2", "FN1"           # near_mTEC_lo (pooled)
)

missing_genes <- setdiff(genes_pooled_check, rownames(samples))
if (length(missing_genes) > 0) {
  warning("Not in object, dropping: ", paste(missing_genes, collapse = ", "))
  genes_pooled_check <- setdiff(genes_pooled_check, missing_genes)
}
stopifnot(length(genes_pooled_check) > 0)

# sanity check: confirm cell_type_unintegrated levels before plotting
cat("Current cell_type_unintegrated levels:\n")
print(table(samples$cell_type_unintegrated))

p_pooled_whole_dataset <- VlnPlot(
  samples, features = genes_pooled_check,
  group.by = "cell_type_unintegrated", pt.size = 0, ncol = 4
) &
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.title.x = element_blank())

print(p_pooled_whole_dataset)
ggsave(paste0(DATE_PREFIX, "-c7-subregion-genes-POOLED-whole-dataset-violin.pdf"),
       p_pooled_whole_dataset, width = 16, height = 16)
ggsave(paste0(DATE_PREFIX, "-c7-subregion-genes-POOLED-whole-dataset-violin.png"),
       p_pooled_whole_dataset, width = 16, height = 16, dpi = 150)

# -----------------------------------------------------------------------------
# Heatmap: top 5 genes from the POOLED vs-OTHER3-internal analysis for each
# c7 sub-region, shown across all whole-dataset cell types -- same style as
# the established Select Markers / PolyKRT / Bautista / Bipolar TEC heatmaps
# -----------------------------------------------------------------------------
genes_pooled_check <- c(
  "CD4", "SPN", "CD96", "RHOH", "CD247",             # near_Nurse (pooled)
  "TNFRSF17", "HAO2", "CRYAB", "SNX22", "CAV1",      # near_cTEC (pooled)
  "MYL9", "CXCL10", "C3", "SLC34A2", "FN1",          # near_mTEC_lo (pooled)
  "FEZF2", "DLX5", "CLDN3", "SRGN", "CXCL5"         # near_mTEC_hi (pooled)
)

missing_genes <- setdiff(genes_pooled_check, rownames(samples))
if (length(missing_genes) > 0) {
  warning("Not in object, dropping: ", paste(missing_genes, collapse = ", "))
  genes_pooled_check <- setdiff(genes_pooled_check, missing_genes)
}
stopifnot(length(genes_pooled_check) > 0)

celltype_order <- c("Nurse", "cTEC", "BMC TEC", "mTEC lo", "mTEC hi",
                    "Mimetic", "Bipolar TEC", "Endothelial", "Pericyte", "Erythroid")

# sanity check: confirm celltype_order matches the object's actual current levels
current_levels <- levels(droplevels(samples$cell_type_unintegrated))
cat("Current cell_type_unintegrated levels:\n")
print(current_levels)
stopifnot(setequal(celltype_order, current_levels))

avg_mat_pooled <- AverageExpression(
  samples, features = genes_pooled_check, assay = "RNA",
  group.by = "cell_type_unintegrated", verbose = TRUE
)$RNA
avg_mat_pooled <- avg_mat_pooled[genes_pooled_check, celltype_order, drop = FALSE]

log_mat_pooled <- log1p(avg_mat_pooled)
z_mat_pooled   <- t(scale(t(log_mat_pooled)))
z_mat_pooled[is.nan(z_mat_pooled)] <- 0

# sanity check: confirm heatmap matrix dimensions match gene count x cell type count
stopifnot(nrow(z_mat_pooled) == length(genes_pooled_check))
stopifnot(ncol(z_mat_pooled) == length(celltype_order))
stopifnot(!any(is.na(z_mat_pooled)))

df_pooled_heatmap <- as.data.frame(z_mat_pooled) %>%
  tibble::rownames_to_column("gene") %>%
  tidyr::pivot_longer(-gene, names_to = "cluster", values_to = "z") %>%
  mutate(
    gene    = factor(gene, levels = rev(genes_pooled_check)),
    cluster = factor(cluster, levels = celltype_order)
  )

p_pooled_heatmap <- ggplot(df_pooled_heatmap, aes(x = cluster, y = gene, fill = z)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient2(
    low = "#F7E27A", mid = "grey95", high = "#1B3B6F",
    midpoint = 0, name = "Relative\nexpr. (z)"
  ) +
  labs(title = "c7 sub-region pooled top-5 genes, across all cell types") +
  theme_minimal(base_size = 10) +
  theme(
    axis.title  = element_blank(),
    plot.title  = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid  = element_blank()
  )

print(p_pooled_heatmap)
ggsave(paste0(DATE_PREFIX, "-c7-subregion-genes-POOLED-heatmap.pdf"),
       p_pooled_heatmap, width = 7, height = 7)

# -----------------------------------------------------------------------------
# Correlation analysis: average expression profile of each c7 sub-region vs.
# every real cell type in the All -- ordered to match the
# established project convention (Nurse, cTEC, BMC TEC, mTEC lo, mTEC hi,
# Mimetic, Endothelial, Bipolar TEC, Pericyte, Erythroid)
# -----------------------------------------------------------------------------
Idents(samples) <- "ident_for_c7_analysis"
# sanity check: confirm Idents matches the c7 sub-region metadata column
stopifnot(all(as.character(Idents(samples)) == unname(samples$ident_for_c7_analysis)))

subregions <- c("near_mTEC_hi", "near_mTEC_lo", "near_cTEC", "near_Nurse")

celltype_order <- c("Nurse", "cTEC", "BMC TEC", "mTEC lo", "mTEC hi",
                    "Mimetic", "Bipolar TEC","Endothelial", "Pericyte", "Erythroid")

# sanity check: confirm celltype_order matches the object's actual current levels
current_levels <- levels(droplevels(samples$cell_type_unintegrated))
cat("Current cell_type_unintegrated levels:\n")
print(current_levels)
stopifnot(setequal(celltype_order, current_levels))

# Average expression, ALL groups in ident_for_c7_analysis
avg_all_groups <- AverageExpression(
  samples, features = VariableFeatures(samples),
  group.by = "ident_for_c7_analysis", assay = "RNA", verbose = TRUE
)$RNA

# AverageExpression() sanitizes "_" to "-" in column names
subregions_sanitized <- c("near-mTEC-hi","near-mTEC-lo","near-cTEC","near-Nurse")
stopifnot(all(subregions_sanitized %in% colnames(avg_all_groups)))

avg_subregions <- avg_all_groups[, subregions_sanitized]
colnames(avg_subregions) <- subregions

# sanity check: confirm restored column names match expected sub-region labels
stopifnot(identical(colnames(avg_subregions), subregions))

# Average expression per real whole-dataset cell type
avg_celltypes <- AverageExpression(
  samples, features = VariableFeatures(samples),
  group.by = "cell_type_unintegrated", assay = "RNA", verbose = TRUE
)$RNA
avg_celltypes <- avg_celltypes[, celltype_order]

# sanity check: confirm same gene set (rows) used for both matrices before correlating
stopifnot(identical(rownames(avg_subregions), rownames(avg_celltypes)))

log_subregions <- as.matrix(log1p(avg_subregions))
log_celltypes  <- as.matrix(log1p(avg_celltypes))

# sanity check: confirm both are now plain numeric matrices
stopifnot(is.matrix(log_subregions), is.numeric(log_subregions))
stopifnot(is.matrix(log_celltypes), is.numeric(log_celltypes))

cor_mat <- cor(log_subregions, log_celltypes, method = "pearson")

# sanity check: confirm correlation matrix dimensions match sub-region x cell-type counts
stopifnot(nrow(cor_mat) == length(subregions))
stopifnot(ncol(cor_mat) == length(celltype_order))

cat("Correlation of each c7 sub-region's average expression profile\n")
cat("against every whole-dataset cell type:\n\n")
print(round(cor_mat, 3))

df_cor <- as.data.frame(cor_mat) %>%
  tibble::rownames_to_column("subregion") %>%
  tidyr::pivot_longer(-subregion, names_to = "cell_type", values_to = "r") %>%
  mutate(
    subregion = factor(subregion, levels = subregions),
    cell_type = factor(cell_type, levels = celltype_order)
  )

p_cor_heatmap <- ggplot(df_cor, aes(x = cell_type, y = subregion, fill = r)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = round(r, 2)), size = 3) +
  scale_fill_gradient2(low = "#F7E27A", mid = "grey95", high = "#1B3B6F",
                       midpoint = median(df_cor$r), name = "Pearson\nr") +
  labs(title = "Bipolar TEC sub-region vs. whole-dataset cell type correlation",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())

print(p_cor_heatmap)
ggsave(paste0(DATE_PREFIX, "-c7-subregion-celltype-correlation-heatmap.pdf"),
       p_cor_heatmap, width = 8, height = 5)




 
























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

