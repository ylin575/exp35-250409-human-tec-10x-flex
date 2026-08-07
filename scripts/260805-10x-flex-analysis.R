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
