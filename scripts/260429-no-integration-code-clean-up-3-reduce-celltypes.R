# exp35-human-tec-10x-flex

# Load Seurat
library(dplyr)
library(Seurat)
library(scCustomize)


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


# add mitochondrial gene reads percentage into metadata
ht67_cd205neg[["percent.mt"]] <- PercentageFeatureSet(ht67_cd205neg, pattern = "^MT-")
ht67_cd205pos[["percent.mt"]] <- PercentageFeatureSet(ht67_cd205pos, pattern = "^MT-")
ht70[["percent.mt"]] <- PercentageFeatureSet(ht70, pattern = "^MT-")
ht71[["percent.mt"]] <- PercentageFeatureSet(ht71, pattern = "^MT-")



# subset and calculate mean and 3*SD for percent.mt (ignore lower bound since
# they tend to be less than 0) and nFeature_RNA

# new method
# meta.data <- bind_rows(ht67_cd205neg@meta.data,
#                        ht67_cd205pos@meta.data,
#                        ht70@meta.data,
#                        ht71@meta.data)
# 
# meta.data %>%
#   group_by(sample.id) %>%
#   summarise(mean_pct_mt = mean(percent.mt),
#             sd_pct_mt = sd(percent.mt),
#             upper.mt = mean_pct_mt + 3 * sd_pct_mt)
#   
# process_seu <- function(seur, min.feature=500, sd.thr=3) {
#   stats <- seur@meta.data %>%
#     summarise(mean_pct_mt = mean(percent.mt),
#               sd_pct_mt = sd(percent.mt),
#               upper.mt = mean_pct_mt + sd.thr * sd_pct_mt,
#               upper.feature=mean(nFeature_RNA) + sd.thr *sd(nFeature_RNA)
#     )
#   subseur <-  subset(seur, subset = percent.mt <= stats$upper.mt)
#   subseur <-  subset(subseur, subset = nFeature_RNA >= min.feature & 
#                        nFeature_RNA <= stats$upper.feature)
#   return(list(stats=stats, seur=subseur))
# }
# 
# ht67_cd205pos_sub <- process_seu(ht67_cd205pos)
# ht67_cd205neg_sub <- process_seu(ht67_cd205neg)


# old method
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

# merge all sample qc stats into one object
#samples <- merge(x = ht67_cd205neg_sub$seur, y = list(ht67_cd205pos_sub, ht70_sub, ht71_sub))

# join layers
samples <- JoinLayers(samples)

# standard analysis steps
samples <- NormalizeData(samples)
samples <- FindVariableFeatures(samples)
samples <- ScaleData(samples)
samples <- RunPCA(samples)


# Examine and visualize PCA results a few different ways
print(samples[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(samples, dims = 1:2, reduction = "pca")
DimPlot(samples, reduction = "pca") + NoLegend()
DimHeatmap(samples, dims = 1:30, cells =50, balanced = TRUE)
ElbowPlot(samples, ndims = 60)

# clustering
samples <- FindNeighbors(samples, dims = 1:18)
samples <- FindClusters(samples, resolution = 0.5)
samples <- RunUMAP(samples, dims = 1:18)


# check some features

DimPlot(samples, reduction = "umap", label = TRUE)
DimPlot(samples, reduction = "umap", label = TRUE, group.by = "sample.id")
DimPlot(samples, reduction = "umap", label = TRUE, split.by = "sample.id")
DimPlot(subset(samples, samples@meta.data$seurat_clusters == 8), reduction = "umap", label = TRUE, split.by = "sample.id")
DimPlot(subset(samples, samples@meta.data$seurat_clusters == 8), reduction = "umap", label = TRUE)
VlnPlot(samples, features = c("TBATA","PRSS16"))

label_size <- 4
order <- TRUE

# cTEC
FeaturePlot_scCustom(samples, reduction = "umap", features = "TBATA", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "PRSS16", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "CLEC2L", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
# Transit-amplifying
FeaturePlot_scCustom(samples, reduction = "umap", features = "MKI67", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "TOP2A", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "CDK1", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "CENPF", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
# mTEC I
FeaturePlot_scCustom(samples, reduction = "umap", features = "CCL19", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "KRT15", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "KRT19", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "FN1", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "CLU", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
# mTEC II
FeaturePlot_scCustom(samples, reduction = "umap", features = "AIRE", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "FEZF2", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
# Post-aire
FeaturePlot_scCustom(samples, reduction = "umap", features = "SPINK5", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
# Ionocyte
FeaturePlot_scCustom(samples, reduction = "umap", features = "CFTR", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "CDH12", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
# aaTEC
FeaturePlot_scCustom(samples, reduction = "umap", features = "VIM", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "IGFBP7", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "SPARC", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "SPARCL1", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "COL1A1", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
# Tuft
FeaturePlot_scCustom(samples, reduction = "umap", features = "POU2F3", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
# Transitioning
FeaturePlot_scCustom(samples, reduction = "umap", features = "TESC", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "CLDN3", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "RGS17", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
# Neuroendocrine
FeaturePlot_scCustom(samples, reduction = "umap", features = "HES6", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "NEUROD1", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
# Neuro
FeaturePlot_scCustom(samples, reduction = "umap", features = "ATOH1", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "USH2A", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
# Muscle
FeaturePlot_scCustom(samples, reduction = "umap", features = "MYOG", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "TTN", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "MYH3", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
# COL4A6-positive TEC
FeaturePlot_scCustom(samples, reduction = "umap", features = "COL4A6", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "COL4A5", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "COL8A1", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
# CCL21-positive TEC
FeaturePlot_scCustom(samples, reduction = "umap", features = "CCL21", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
# Nurse
FeaturePlot_scCustom(samples, reduction = "umap", features = "TBATA", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "PRSS16", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "CD3E", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "PTPRC", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
# Ciliated
FeaturePlot_scCustom(samples, reduction = "umap", features = "FOXJ1", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
# CFC1-positive cTECs
FeaturePlot_scCustom(samples, reduction = "umap", features = "TBATA", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "PRSS16", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "CFC1", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "TP53AIP1", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "UBD", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
# All TECs
FeaturePlot_scCustom(samples, reduction = "umap", features = "TP63", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
# Endothelial
FeaturePlot_scCustom(samples, reduction = "umap", features = "PECAM1", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = F)
# Other markers
FeaturePlot_scCustom(samples, reduction = "umap", features = "LY75", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "PRSS16", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "PSMB11", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "COL4A5", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "COL4A6", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "ITGB4", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "CD200", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "ITGA6", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "BCAM", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "EPCAM", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "CD24", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "CD74", colors_use = viridis_light_high, order = order)


# save object before annotation
saveRDS(samples, file = "data/rds-objects/260428-before-annotation.rds")
samples <- readRDS(file = "data/rds-objects/260428-before-annotation.rds")


# annotate cell types
test <- samples
new.cluster.ids <- c("cTEC",
                     "cTEC",
                     "mTEC I",
                     "Nurse",
                     "COL4A6+ TEC",
                     "COL4A6+ TEC",
                     "mTEC I",
                     "mTEC II",
                     "mTEC I",
                     "mTEC I",
                     "COL4A6+ TEC",
                     "TD-TEC",
                     "TD-TEC",
                     "TD-TEC",
                     "Nurse",
                     "TD-TEC",
                     "Nurse",
                     "mTEC I",
                     "mTEC II",
                     "TD-TEC",
                     "TD-TEC",
                     "Endothelial",
                     "cTEC",
                     "cTEC",
                     "mTEC I",
                     "cTEC",
                     "mTEC I")
levels(test)
names(new.cluster.ids) <- levels(test)
names(new.cluster.ids)
new.cluster.ids
View(new.cluster.ids)
test <- RenameIdents(test, new.cluster.ids)
DimPlot(test, reduction = "umap", label = FALSE, pt.size = 0.4)
test@meta.data$cell_type <- Idents(test)
DimPlot(test, reduction = "umap", label = FALSE, group.by = "cell_type", split.by = "sample.id")
#new colors for cell clusters: palman values
DimPlot(test, reduction = "umap", label = FALSE,
        cols = c("#D8A767","#F47D2B","#D24B27","#3D3D3D","#7E1416","#0570B0",
                 "#89288F","#3BBCA8","#B45F06","#89C75F","#FEE500","#208A42",
                 "#f30892")
)


# save object after annotation
saveRDS(test, file = "data/rds-objects/260429-after-annotation-1.rds")
test <- readRDS("data/rds-objects/260429-after-annotation-1.rds")



