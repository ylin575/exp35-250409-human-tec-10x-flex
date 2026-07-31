# exp35-human-tec-10x-flex

# Load Seurat
library(dplyr)
library(Seurat)
library(scCustomize)
library(ggplot2)

DATE_PREFIX <- format(Sys.time(), "%y%m%d-%H%M")

# load the count matrices
ht67_cd205neg.data <- Read10X_h5("rawdata/h5-files/ht67-cd205neg/sample_filtered_feature_bc_matrix.h5", 
                                 use.names = TRUE, unique.features = TRUE)
ht67_cd205pos.data <- Read10X_h5("rawdata/h5-files/ht67-cd205pos/sample_filtered_feature_bc_matrix.h5", 
                                 use.names = TRUE, unique.features = TRUE)
ht70.data <- Read10X_h5("rawdata/h5-files/ht70/sample_filtered_feature_bc_matrix.h5", 
                                 use.names = TRUE, unique.features = TRUE)
ht71.data <- Read10X_h5("rawdata/h5-files/ht71/sample_filtered_feature_bc_matrix.h5", 
                                 use.names = TRUE, unique.features = TRUE)

# Initialize the Seurat object with the raw (non-normalized data) counts
ht67_cd205neg <- CreateSeuratObject(counts = ht67_cd205neg.data, project = "htec-10xflex", 
                                    min.cells = 3,
                                    min.features = 200)
ht67_cd205pos <- CreateSeuratObject(counts = ht67_cd205pos.data, project = "htec-10xflex", 
                                    min.cells = 3,
                                    min.features = 200)
ht70 <- CreateSeuratObject(counts = ht70.data, project = "htec-10xflex", 
                                    min.cells = 3,
                                    min.features = 200)
ht71 <- CreateSeuratObject(counts = ht71.data, project = "htec-10xflex", 
                                    min.cells = 3,
                                    min.features = 200)

# add sample ID to metadata
ht67_cd205neg[["sample.id"]] <- "ht67-cd205neg"
ht67_cd205pos[["sample.id"]] <- "ht67-cd205pos"
ht70[["sample.id"]] <- "ht70"
ht71[["sample.id"]] <- "ht71"

# add donor to metadata (HT67's two CD205 gates are the same individual)
ht67_cd205neg[["donor"]] <- "HT67"
ht67_cd205pos[["donor"]] <- "HT67"
ht70[["donor"]] <- "HT70"
ht71[["donor"]] <- "HT71"

# add mitochondrial gene reads percentage into metadata
ht67_cd205neg[["percent.mt"]] <- PercentageFeatureSet(ht67_cd205neg, pattern = "^MT-")
ht67_cd205pos[["percent.mt"]] <- PercentageFeatureSet(ht67_cd205pos, pattern = "^MT-")
ht70[["percent.mt"]] <- PercentageFeatureSet(ht70, pattern = "^MT-")
ht71[["percent.mt"]] <- PercentageFeatureSet(ht71, pattern = "^MT-")



# QC filter
ht67.cd205neg.mean.percent.mt <- mean(ht67_cd205neg@meta.data$percent.mt)
ht67.cd205neg.sd.percent.mt <- sd(ht67_cd205neg@meta.data$percent.mt)
ht67.cd205neg.upper.mt <- ht67.cd205neg.mean.percent.mt + 3 * ht67.cd205neg.sd.percent.mt
ht67.cd205neg.upper.feature <- mean(ht67_cd205neg@meta.data$nFeature_RNA) + 3*sd(ht67_cd205neg@meta.data$nFeature_RNA)
ht67_cd205neg_sub <- subset(ht67_cd205neg, subset = percent.mt <= ht67.cd205neg.upper.mt)
ht67_cd205neg_sub <- subset(ht67_cd205neg_sub, subset = nFeature_RNA >= 500 & nFeature_RNA <= ht67.cd205neg.upper.feature)


ht67.cd205pos.mean.percent.mt <- mean(ht67_cd205pos@meta.data$percent.mt)
ht67.cd205pos.sd.percent.mt <- sd(ht67_cd205pos@meta.data$percent.mt)
ht67.cd205pos.upper.mt <- ht67.cd205pos.mean.percent.mt + 3 * ht67.cd205pos.sd.percent.mt
ht67.cd205pos.upper.feature <- mean(ht67_cd205pos@meta.data$nFeature_RNA) + 3*sd(ht67_cd205pos@meta.data$nFeature_RNA)
ht67_cd205pos_sub <- subset(ht67_cd205pos, subset = percent.mt <= ht67.cd205pos.upper.mt)
ht67_cd205pos_sub <- subset(ht67_cd205pos_sub, subset = nFeature_RNA >= 500 & nFeature_RNA <= ht67.cd205pos.upper.feature)


ht70.mean.percent.mt <- mean(ht70@meta.data$percent.mt)
ht70.sd.percent.mt <- sd(ht70@meta.data$percent.mt)
ht70.upper.mt <- ht70.mean.percent.mt + 3 * ht70.sd.percent.mt
ht70.upper.feature <- mean(ht70@meta.data$nFeature_RNA) + 3*sd(ht70@meta.data$nFeature_RNA)
ht70_sub <- subset(ht70, subset = percent.mt <= ht70.upper.mt)
ht70_sub <- subset(ht70_sub, subset = nFeature_RNA >= 500 & nFeature_RNA <= ht70.upper.feature)


ht71.mean.percent.mt <- mean(ht71@meta.data$percent.mt)
ht71.sd.percent.mt <- sd(ht71@meta.data$percent.mt)
ht71.upper.mt <- ht71.mean.percent.mt + 3 * ht71.sd.percent.mt
ht71.upper.feature <- mean(ht71@meta.data$nFeature_RNA) + 3*sd(ht71@meta.data$nFeature_RNA)
ht71_sub <- subset(ht71, subset = percent.mt <= ht71.upper.mt)
ht71_sub <- subset(ht71_sub, subset = nFeature_RNA >= 500 & nFeature_RNA <= ht71.upper.feature)


# merge all samples into one object
samples <- merge(x = ht67_cd205neg_sub, y = list(ht67_cd205pos_sub, ht70_sub, 
                                                 ht71_sub))
# join layers
samples <- JoinLayers(samples)

# sanity check: verify donor grouping with sample.id
table(samples$sample.id, samples$donor)

# standard analysis steps (moved BEFORE Harmony -- PCA must exist first)
samples <- NormalizeData(samples)
samples <- FindVariableFeatures(samples)
samples <- ScaleData(samples)
samples <- RunPCA(samples)

# Examine and visualize PCA results a few different ways
print(samples[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(samples, dims = 1:2, reduction = "pca")
DimPlot(samples, reduction = "pca") + NoLegend()
DimHeatmap(samples, dims = 1:30, cells = 50, balanced = TRUE)
ElbowPlot(samples, ndims = 60)

# clustering
samples <- FindNeighbors(samples, dims = 1:18)
samples <- FindClusters(samples, resolution = 0.5)
samples <- RunUMAP(samples, dims = 1:18)

# Cluster count and size distribution
cat("Clusters found:", length(unique(samples$seurat_clusters)), "\n")
print(table(samples$seurat_clusters))

# DimPlots
DimPlot(samples, reduction = "umap", label = TRUE)
DimPlot(samples, reduction = "umap", label = TRUE, group.by = "sample.id")
DimPlot(samples, reduction = "umap", label = TRUE, split.by = "sample.id")
DimPlot(subset(samples, samples@meta.data$seurat_clusters == 6), 
        reduction = "umap", label = TRUE)
DimPlot(subset(samples, samples@meta.data$seurat_clusters == 6), 
        reduction = "umap", label = TRUE, split.by = "sample.id")
DimPlot_scCustom(samples, cells.highlight = WhichCells(samples, idents = "2"))



save_features <- function(feature, seur, reduction = "umap", label = TRUE, 
                          label_size = 4, order = TRUE, split.by = NULL) {
  
  p <- FeaturePlot_scCustom(seur, reduction = reduction, features = feature,
                            colors_use = viridis_light_high, label = label, 
                            label.size = label_size, order = order, 
                            split.by = split.by)
  
  ggsave(plot=p, paste0("results/",DATE_PREFIX,"-featureplot-",obj_name,feature, ".pdf"),
         width=8, height=7)
}


# batch save featureplots of interest
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

save_features <- function(feature, seur, reduction = "umap", label = TRUE, 
                          label_size = 4, order = TRUE, split.by = NULL) {
  
  p <- FeaturePlot_scCustom(seur, reduction = reduction, features = feature,
                            colors_use = viridis_light_high, label = label, 
                            label.size = label_size, order = order, 
                            split.by = split.by)
  
  ggsave(plot=p, paste0("results/",DATE_PREFIX,"-featureplot-",obj_name,feature, ".pdf"),
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


# Other markers
FeaturePlot_scCustom(samples, reduction = "umap", features = "CD200", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "ITGA6", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "BCAM", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "EPCAM", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "CD74", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "CD200", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "FOXN1", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "CCL21", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "PDPN", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "DLK2", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "CCL19", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "KRT13", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "KRT15", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "FN1", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "COL6A3", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "ITGA6", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "ITGA5", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "IFITM3", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "LIFR", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "CD34", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "CCL2", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "CCN2", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "IGFBP6", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "UBD", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "ITGAX", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "FAP", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "COL1A1", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "COL1A2", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "COL2A1", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "CCL25", colors_use = viridis_light_high, order = order)

# cytokeratin genes
krt <- c("KRT1","KRT2","KRT3","KRT4","KRT5","KRT6A","KRT6B","KRT6C","KRT7",
         "KRT8","KRT9","KRT10","KRT12","KRT13","KRT14","KRT15","KRT16","KRT17",
         "KRT19","KRT20","KRT23","KRT24","KRT25","KRT27","KRT28","KRT30",
         "KRT31","KRT32","KRT33A","KRT33B","KRT34","KRT35","KRT36","KRT37",
         "KRT38","KRT39","KRT40","KRT71","KRT72","KRT73","KRT74","KRT75",
         "KRT76","KRT77","KRT78","KRT79","KRT80")


VlnPlot(samples, features = c("PSMB11", "PTPRC"), idents = c("4","13","14","15","19","21"))
VlnPlot(samples, features = c("CD4", "CD8A"), idents = c("4", "13", "14", "15"))


# save object before annotation
saveRDS(samples, file = paste0("data/rds-objects/",DATE_PREFIX,"-before-annotation.rds"))



##############################################
# Find Markers #
##############################################
DATE_PREFIX <- format(Sys.time(), "%y%m%d-%H%M")

Idents(samples) <- "seurat_clusters"

all_markers <- FindAllMarkers(
  samples,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25,
  verbose         = TRUE
)

# Sanity check: every cluster represented, no unexpected gaps
sanity_clusters_found  <- sort(unique(all_markers$cluster))
sanity_clusters_expect <- sort(unique(samples$seurat_clusters))
cat("Clusters with markers found:", length(sanity_clusters_found), "\n")
cat("Clusters expected:", length(sanity_clusters_expect), "\n")
missing_clusters <- setdiff(sanity_clusters_expect, sanity_clusters_found)
if (length(missing_clusters) > 0) {
  cat("WARNING: no markers passed threshold for cluster(s):",
      paste(missing_clusters, collapse = ", "), "\n")
} else {
  cat("PASS: all clusters have at least one marker gene.\n")
}

# Top 15 per cluster, for a quick scan
top_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 15) %>%
  ungroup()

write.csv(all_markers, paste0(DATE_PREFIX, "-unintegrated-all-markers-FULL.csv"),
          row.names = FALSE)
write.csv(top_markers, paste0(DATE_PREFIX, "-unintegrated-all-markers-TOP15.csv"),
          row.names = FALSE)


# Cluster 6 specifically, since that's the population of interest
cat("\n=== Cluster 6 top 15 markers ===\n")
print(top_markers %>% filter(cluster == "6") %>%
        select(gene, avg_log2FC, pct.1, pct.2, p_val_adj))


#######################################################


# assess the former cluster 8 that was previously spread widely

# cTEC/mTEC-I marker panels (same as used throughout the C1-C6 analysis)
cTEC_markers   <- c("PSMB11", "PRSS16", "LY75")
mTEC_I_markers <- c("KRT14", "CCL19", "KRT15", "KRT19")

cTEC_markers   <- intersect(cTEC_markers, rownames(samples))
mTEC_I_markers <- intersect(mTEC_I_markers, rownames(samples))

# Strict double-positive test, same detected-gene-count method used before
counts_now <- GetAssayData(samples, assay = "RNA", layer = "counts")
cTEC_detected   <- Matrix::colSums(counts_now[cTEC_markers, , drop = FALSE] > 0)
mTEC_I_detected <- Matrix::colSums(counts_now[mTEC_I_markers, , drop = FALSE] > 0)

samples$double_pos_strict <- (cTEC_detected   >= ceiling(length(cTEC_markers)   / 2)) &
  (mTEC_I_detected >= ceiling(length(mTEC_I_markers) / 2))

# Per-cluster summary, ranked
cluster_check <- samples@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(
    n = n(),
    pct_double_pos = round(100 * mean(double_pos_strict), 1),
    .groups = "drop"
  ) %>%
  arrange(desc(pct_double_pos))

print(cluster_check)



#######################################################



# annotate post harmony integrated cell types
test <- samples
new.cluster.ids <- c("cTEC",
                     "mTEC lo",
                     "cTEC",
                     "BMC TEC",
                     "Nurse",
                     "mTEC lo",
                     "BP TEC",
                     "mTEC hi",
                     "Mimetic",
                     "mTEC lo",
                     "Mimetic",
                     "Mimetic",
                     "Mimetic",
                     "Nurse",
                     "Nurse",
                     "Nurse",
                     "Pericyte",
                     "Mimetic",
                     "Mimetic",
                     "Erythroid",
                     "Mimetic",
                     "Endothelial")
levels(test)
names(new.cluster.ids) <- levels(test)
names(new.cluster.ids)
new.cluster.ids
test <- RenameIdents(test, new.cluster.ids)
DimPlot(test, reduction = "umap", label = FALSE, pt.size = 0.4)
DimPlot(test, reduction = "umap", label = TRUE, pt.size = 0.4)
test@meta.data$cell_type <- Idents(test)
DimPlot(test, reduction = "umap", label = FALSE, group.by = "cell_type", split.by = "sample.id")
#new colors for cell clusters: palman values
DimPlot(test, reduction = "umap", label = FALSE,
        cols = c("#D8A767","#F47D2B","#D24B27","#E7298A","#7E1416","#0570B0",
                 "#89C75F","#3D3D3D","#208A42","#D51F26")
)


# save object after annotation
saveRDS(test, file = paste0("data/rds-objects/",DATE_PREFIX,"-after-annotation.rds"))
#tec <- readRDS("data/rds-objects/260501-after-annotation-1.rds")


# 260706 plot for revision JI
genes <- c("PTPRC","PSMB11","PRSS16","LY75","COL4A5","COL4A6","ITGB4","CD200",
           "CLU","AIRE","FEZF2","CD24","PECAM1")

tec2 <- tec
my_order <- c("Nurse","cTEC","BMC+ TEC","mTEC I","mTEC II","Mimetic","Endothelial")
Idents(tec2) <- factor(Idents(tec2), levels = my_order)
VlnPlot(tec2, ncol = 3, features = genes)



# assess cell cycle score for the nurse cell subsets
samples <- CellCycleScoring(
  samples,
  s.features   = cc.genes.updated.2019$s.genes,
  g2m.features = cc.genes.updated.2019$g2m.genes,
  set.ident    = FALSE
)

cat("=== Cell cycle phase by nurse subcluster ===\n")
print(table(samples$Phase, samples$seurat_clusters)[, c("4","13","14","15")])



# assess potential doublet contamination

library(scDblFinder)
library(SingleCellExperiment)

# scDblFinder must respect GEM-well boundaries -- doublets only form within
# a single sample, not across samples. Use sample.id (not donor) here, since
# that's the actual physical unit of capture.
sce <- as.SingleCellExperiment(samples, assay = "RNA")
set.seed(42)
sce <- scDblFinder(sce, samples = "sample.id", verbose = TRUE)

samples$scDblFinder_class <- as.character(sce$scDblFinder.class)
samples$scDblFinder_score <- as.numeric(sce$scDblFinder.score)

# Sanity check
stopifnot(!any(is.na(samples$scDblFinder_class)))
cat("Overall doublet rate:", round(100 * mean(samples$scDblFinder_class == "doublet"), 1), "%\n")

# Per-cluster doublet rate, all clusters for context, 19 and 17 highlighted
library(dplyr)
doublet_by_cluster <- samples@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(
    n = n(),
    pct_doublet = round(100 * mean(scDblFinder_class == "doublet"), 1),
    .groups = "drop"
  ) %>%
  arrange(desc(pct_doublet))

print(doublet_by_cluster)

cat("\n=== Clusters 17 and 19 specifically ===\n")
print(doublet_by_cluster %>% filter(seurat_clusters %in% c("17", "19")))

