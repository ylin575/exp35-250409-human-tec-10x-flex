# exp35-human-tec-10x-flex

# Load Packages
library(dplyr)
library(Seurat)
library(patchwork)
library(scIntegrationMetrics)


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


# > length(rownames(ht67_cd205neg.data))
# [1] 18082
# > length(rownames(ht67_cd205neg))
# [1] 17237
# > length(rownames(ht67_cd205pos.data))
# [1] 18082
# > length(rownames(ht67_cd205pos))
# [1] 15816
# > length(rownames(ht70.data))
# [1] 18082
# > length(rownames(ht70))
# [1] 17007
# > length(rownames(ht71.data))
# [1] 18082
# > length(rownames(ht71))
# [1] 16958


saveRDS(ht67_cd205neg, file = "data/rds-objects/250414-ht67-cd205neg-before-QC.rds")
saveRDS(ht67_cd205pos, file = "data/rds-objects/250414-ht67-cd205pos-before-QC.rds")
saveRDS(ht70, file = "data/rds-objects/250414-ht70-before-QC.rds")
saveRDS(ht71, file = "data/rds-objects/250414-ht71-before-QC.rds")

ht67_cd205neg <- readRDS(file = "data/rds-objects/250414-ht67-cd205neg-before-QC.rds")
ht67_cd205pos <- readRDS(file = "data/rds-objects/250414-ht67-cd205pos-before-QC.rds")
ht70 <- readRDS(file = "data/rds-objects/250414-ht70-before-QC.rds")
ht71 <- readRDS(file = "data/rds-objects/250414-ht71-before-QC.rds")

# add sample ID to metadata
ht67_cd205neg[["sample.id"]] <- "ht67-cd205neg"
ht67_cd205pos[["sample.id"]] <- "ht67-cd205pos"
ht70[["sample.id"]] <- "ht70"
ht71[["sample.id"]] <- "ht71"

# perform QC

ht67_cd205neg[["percent.mt"]] <- PercentageFeatureSet(ht67_cd205neg, pattern = "^MT-")
ht67_cd205pos[["percent.mt"]] <- PercentageFeatureSet(ht67_cd205pos, pattern = "^MT-")
ht70[["percent.mt"]] <- PercentageFeatureSet(ht70, pattern = "^MT-")
ht71[["percent.mt"]] <- PercentageFeatureSet(ht71, pattern = "^MT-")



# original dataset dimensions
dim(ht67_cd205neg)
dim(ht67_cd205pos)
dim(ht70)
dim(ht71)

# subset and calculate mean and 3*SD for percent.mt (ignore lower bound since
# they tend to be less than 0) and nFeature_RNA


meta.data <- bind_rows(ht67_cd205neg@meta.data,
                       ht67_cd205pos@meta.data,
                       ht70@meta.data,
                       ht71@meta.data)

meta.data %>%
  group_by(sample.id) %>%
  summarise(mean_pct_mt = mean(percent.mt),
            sd_pct_mt = sd(percent.mt),
            upper.mt = mean_pct_mt + 3 * sd_pct_mt)
  
process_seu <- function(seur, min.feature=500, sd.thr=3) {
  stats <- seur@meta.data %>%
    summarise(mean_pct_mt = mean(percent.mt),
              sd_pct_mt = sd(percent.mt),
              upper.mt = mean_pct_mt + sd.thr * sd_pct_mt,
              upper.feature=mean(nFeature_RNA) + sd.thr *sd(nFeature_RNA)
    )
  subseur <-  subset(seur, subset = percent.mt <= stats$upper.mt)
  subseur <-  subset(subseur, subset = nFeature_RNA >= min.feature & 
                       nFeature_RNA <= stats$upper.feature)
  return(list(stats=stats, seur=subseur))
}

ht67_cd205pos_sub <- process_seu(ht67_cd205pos)
ht67_cd205neg_sub <- process_seu(ht67_cd205neg)



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


# Visualize QC metrics as a violin plot
VlnPlot(ht67_cd205neg_sub, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(ht67_cd205pos_sub, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(ht70_sub, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(ht71_sub, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(ht67_cd205neg, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ht67_cd205neg, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(ht67_cd205pos, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ht67_cd205pos, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(ht70, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ht70, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(ht71, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ht71, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


plot1 <- FeatureScatter(ht67_cd205neg_sub, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ht67_cd205neg_sub, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(ht67_cd205pos_sub, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ht67_cd205pos_sub, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(ht70_sub, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ht70_sub, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(ht71_sub, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ht71_sub, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


###########################
# mean and median pre qc filter

mean(ht67_cd205neg@meta.data$nFeature_RNA)
mean(ht67_cd205pos@meta.data$nFeature_RNA)
mean(ht70@meta.data$nFeature_RNA)
mean(ht71@meta.data$nFeature_RNA)

mean(ht67_cd205neg@meta.data$nCount_RNA)
mean(ht67_cd205pos@meta.data$nCount_RNA)
mean(ht70@meta.data$nCount_RNA)
mean(ht71@meta.data$nCount_RNA)

median(ht67_cd205neg_sub@meta.data$nFeature_RNA)
median(ht67_cd205pos_sub@meta.data$nFeature_RNA)
median(ht70_sub@meta.data$nFeature_RNA)
median(ht71_sub@meta.data$nFeature_RNA)

median(ht67_cd205neg_sub@meta.data$nCount_RNA)
median(ht67_cd205pos_sub@meta.data$nCount_RNA)
median(ht70_sub@meta.data$nCount_RNA)
median(ht71_sub@meta.data$nCount_RNA)

###########################
# mean and median post qc filter

mean(ht67_cd205neg_sub@meta.data$nFeature_RNA)
mean(ht67_cd205pos_sub@meta.data$nFeature_RNA)
mean(ht70_sub@meta.data$nFeature_RNA)
mean(ht71_sub@meta.data$nFeature_RNA)

mean(ht67_cd205neg_sub@meta.data$nCount_RNA)
mean(ht67_cd205pos_sub@meta.data$nCount_RNA)
mean(ht70_sub@meta.data$nCount_RNA)
mean(ht71_sub@meta.data$nCount_RNA)

median(ht67_cd205neg_sub@meta.data$nFeature_RNA)
median(ht67_cd205pos_sub@meta.data$nFeature_RNA)
median(ht70_sub@meta.data$nFeature_RNA)
median(ht71_sub@meta.data$nFeature_RNA)

median(ht67_cd205neg_sub@meta.data$nCount_RNA)
median(ht67_cd205pos_sub@meta.data$nCount_RNA)
median(ht70_sub@meta.data$nCount_RNA)
median(ht71_sub@meta.data$nCount_RNA)





# # check original dimensions of the subsetted objects
dim(ht67_cd205neg)
dim(ht67_cd205pos)
dim(ht70)
dim(ht71)
dim(ht67_cd205neg_sub)
dim(ht67_cd205pos_sub)
dim(ht70_sub)
dim(ht71_sub)
# > dim(ht67_cd205neg)
# [1] 17237  8280
# > dim(ht67_cd205pos)
# [1] 15816 12953
# > dim(ht70)
# [1] 17007 13246
# > dim(ht71)
# [1] 16958 11073
# > dim(ht67_cd205neg_sub)
# [1] 17237  8158
# > dim(ht67_cd205pos_sub)
# [1] 15816 12662
# > dim(ht70_sub)
# [1] 17007 13043
# > dim(ht71_sub)
# [1] 16958 10703


# merge all samples into one object
samples <- merge(x = ht67_cd205neg_sub, y = list(ht67_cd205pos_sub, ht70_sub, 
                                               ht71_sub))

samples <- merge(x = ht67_cd205neg_sub$seur, y = list(ht67_cd205pos_sub, ht70_sub, 
                                                 ht71_sub))
# standard analysis steps
samples <- NormalizeData(samples)
samples <- FindVariableFeatures(samples)
samples <- ScaleData(samples)
samples <- RunPCA(samples, features = VariableFeatures(object = samples))

# Examine and visualize PCA results a few different ways
print(samples[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(samples, dims = 1:2, reduction = "pca")
DimPlot(samples, reduction = "pca") + NoLegend()
DimHeatmap(samples, dims = 1:15, cells = 6, balanced = TRUE)
ElbowPlot(samples)

# clustering
samples <- FindNeighbors(samples, dims = 1:18)
samples <- FindClusters(samples, resolution = 0.9, 
                      cluster.name = "unintegrated_clusters")

# run non-linear dimensional reduction
samples <- RunUMAP(samples, dims = 1:18, 
                 reduction.name = "umap.unintegrated")

# plot umap before integration
DimPlot(samples, reduction = "umap.unintegrated")
# group by sample id
DimPlot(samples, reduction = "umap.unintegrated", group.by = "sample.id")
# split by sample id
DimPlot(samples, reduction = "umap.unintegrated", split.by = "sample.id")

# check some features
FeaturePlot(samples, reduction = "umap.unintegrated", features = "AIRE", split.by = "sample.id")
FeaturePlot(samples, reduction = "umap.unintegrated", features = "PRSS16", split.by = "sample.id")
FeaturePlot(samples, reduction = "umap.unintegrated", features = "PSMB11", split.by = "sample.id")
FeaturePlot(samples, reduction = "umap.unintegrated", features = "LY75", split.by = "sample.id")
FeaturePlot(samples, reduction = "umap.unintegrated", features = "COL4A6", split.by = "sample.id")
FeaturePlot(samples, reduction = "umap.unintegrated", features = "ITGA6", split.by = "sample.id")
FeaturePlot(samples, reduction = "umap.unintegrated", features = "THY1", split.by = "sample.id")
FeaturePlot(samples, reduction = "umap.unintegrated", features = "BCAM", split.by = "sample.id")
FeaturePlot(samples, reduction = "umap.unintegrated", features = "CLEC2L", split.by = "sample.id")




#save before integration
saveRDS(samples, file = "data/rds-objects/250415-before-integration.rds")
samples <- readRDS("data/rds-objects/250415-before-integration.rds")


# perform integration on the merged samples; wi = with integration
sampleswi <- IntegrateLayers(object = samples, method = CCAIntegration, 
                           orig.reduction = "pca", 
                           new.reduction = "integrated.cca",
                           verbose = FALSE)
# 250415 took 43m 58s


# save object after integrating, before joining layers
saveRDS(sampleswi, file = "data/rds-objects/250415-after-integration.rds")
sampleswi <- readRDS(file = "data/rds-objects/250415-after-integration.rds")


# rejoin the layers after integration
sampleswi[["RNA"]] <- JoinLayers(sampleswi[["RNA"]])


# remove objects
rm(ht67_cd205neg, ht67_cd205pos, ht70, ht71)
rm(ht67_cd205neg_sub, ht67_cd205pos_sub, ht70_sub, ht71_sub)


# examine and visualize PCA results a few different ways
print(sampleswi[["integrated.cca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(sampleswi, dims = 1:2, reduction = "integrated.cca")
# compare DimPlot from before and after integration
DimPlot(samples, reduction = "pca") + NoLegend()
DimPlot(sampleswi, reduction = "integrated.cca") + NoLegend()
DimHeatmap(sampleswi, dims = 1:15, cells = 6, balanced = TRUE)
ElbowPlot(sampleswi, reduction = "integrated.cca")


# find clusters post integration and joining
sampleswi <- FindNeighbors(sampleswi, reduction = "integrated.cca", dims = 1:18)
sampleswi <- FindClusters(sampleswi, resolution = 0.9)

# run non-linear dimensional reduction
sampleswi <- RunUMAP(sampleswi, dims = 1:18, reduction = "integrated.cca")

# plot graphs before combining sample IDs
DimPlot(sampleswi, reduction = "umap")
DimPlot(sampleswi, reduction = "umap", group.by = "sample.id")
DimPlot(sampleswi, reduction = "umap", split.by = "sample.id")
DimPlot(samples, reduction = "umap.unintegrated", split.by = "sample.id")

# plot unintegrated samples
FeaturePlot(samples, reduction = "umap.unintegrated", features = "AIRE", split.by = "sample.id")
FeaturePlot(samples, reduction = "umap.unintegrated", features = "PRSS16", split.by = "sample.id")
FeaturePlot(samples, reduction = "umap.unintegrated", features = "PSMB11", split.by = "sample.id")
FeaturePlot(samples, reduction = "umap.unintegrated", features = "LY75", split.by = "sample.id")
FeaturePlot(samples, reduction = "umap.unintegrated", features = "COL4A6", split.by = "sample.id")
FeaturePlot(samples, reduction = "umap.unintegrated", features = "TTN", split.by = "sample.id")
FeaturePlot(samples, reduction = "umap.unintegrated", features = "KRT5", split.by = "sample.id")
FeaturePlot(samples, reduction = "umap.unintegrated", features = "KRT8", split.by = "sample.id")
FeaturePlot(samples, reduction = "umap.unintegrated", features = "KRT14", split.by = "sample.id")
FeaturePlot(samples, reduction = "umap.unintegrated", features = "KRT15", split.by = "sample.id")
FeaturePlot(samples, reduction = "umap.unintegrated", features = "ITGA6", split.by = "sample.id")
FeaturePlot(samples, reduction = "umap.unintegrated", features = "BCAM", split.by = "sample.id")
FeaturePlot(samples, reduction = "umap.unintegrated", features = "PTPRC", split.by = "sample.id")
FeaturePlot(samples, reduction = "umap.unintegrated", features = "ITGAX", split.by = "sample.id")
FeaturePlot(samples, reduction = "umap.unintegrated", features = "CD3E", split.by = "sample.id")

# plot integrated samples
FeaturePlot(sampleswi, reduction = "umap", features = "AIRE", split.by = "sample.id")
FeaturePlot(sampleswi, reduction = "umap", features = "PRSS16", split.by = "sample.id")
FeaturePlot(sampleswi, reduction = "umap", features = "PSMB11", split.by = "sample.id")
FeaturePlot(sampleswi, reduction = "umap", features = "LY75", split.by = "sample.id")
FeaturePlot(sampleswi, reduction = "umap", features = "COL4A6", split.by = "sample.id")
FeaturePlot(sampleswi, reduction = "umap", features = "MYOG", split.by = "sample.id")
FeaturePlot(sampleswi, reduction = "umap", features = "KRT5", split.by = "sample.id")
FeaturePlot(sampleswi, reduction = "umap", features = "KRT8", split.by = "sample.id")
FeaturePlot(sampleswi, reduction = "umap", features = "KRT14", split.by = "sample.id")
FeaturePlot(sampleswi, reduction = "umap", features = "KRT15", split.by = "sample.id")


# post-integration qc plots
plot1 <- FeatureScatter(sampleswi, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "sample.id")
plot2 <- FeatureScatter(sampleswi, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "sample.id")
plot1 + plot2


# evaluate integration with iLISI














#############################################################

# codes below not used, from another script, used for reference





ht67_cd205neg@meta.data$orig.ident
Idents(ht67_cd205neg)

ht67_cd205neg
str(ht67_cd205neg)

# lets examine a few gens in the first thirty cells
head(ht67_cd205neg.data)
ht67_cd205neg.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]


# ===========================================================


#> QC metrics
#> 1. the number of unique genes detected in each cell
#>     a. low quality cells/empty droplets often have very few genes
#>     b. cell multiplets may exhibit abnormally high gene count
#> 2. the total number of molecules detected within a cell
#> 3. the percent of reads that map to the mitochondrial genome
#>     a. low quality / dying cells often exhibit extensive mitochondrial
#>        contamination

# The [[ operator can add columns to object metadata. This is a great place to 
# stash QC stats
ht67_cd205neg[["percent.mt"]] <- PercentageFeatureSet(ht67_cd205neg, pattern = "^MT-")

# QC metrics are stored in @meta.data slot
# Show QC metrics for the first 5 cells
head(ht67_cd205neg@meta.data, 5)

# Visualize QC metrics as a violin plot
VlnPlot(ht67_cd205neg, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(ht67_cd205neg, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ht67_cd205neg, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# select only cells that meet the following criteria
ht67_cd205neg <- subset(ht67_cd205neg, subset = nFeature_RNA > 200 & nFeature_RNA < 8750 & percent.mt < 2.5)


# ===========================================================


# data normalization
ht67_cd205neg <- NormalizeData(ht67_cd205neg, normalization.method = "LogNormalize", scale.factor = 10000)


# ===========================================================


# find highly variable genes across cells and use these for downstream procedures
ht67_cd205neg <- FindVariableFeatures(ht67_cd205neg, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(ht67_cd205neg), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(ht67_cd205neg)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


# ===========================================================


# scale data needed prior to running PCA
# Shifts the expression of each gene, so that the mean expression across cells is 0
# Scales the expression of each gene, so that the variance across cells is 1
# This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
# The results of this are stored in ht67_cd205neg[["RNA"]]$scale.data
# By default, only variable features are scaled.
# You can specify the features argument to scale additional features
all.genes <- rownames(ht67_cd205neg)
ht67_cd205neg <- ScaleData(ht67_cd205neg, features = all.genes)


# ===========================================================


# perform linear dimensional reduction
ht67_cd205neg <- RunPCA(ht67_cd205neg, features = VariableFeatures(object = ht67_cd205neg))

# Seurat provides several useful ways of visualizing both cells and features that 
# define the PCA, including VizDimReduction(), DimPlot(), and DimHeatmap()

# Examine and visualize PCA results a few different ways
print(ht67_cd205neg[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(ht67_cd205neg, dims = 1:2, reduction = "pca")
DimPlot(ht67_cd205neg, reduction = "pca") + NoLegend()
DimHeatmap(ht67_cd205neg, dims = 1:5, cells = 6, balanced = TRUE)
ElbowPlot(ht67_cd205neg)


# ===========================================================


# cluster cells
# dims is the number of PCA dimensions to use, based on examining pca results above
ht67_cd205neg <- FindNeighbors(ht67_cd205neg, dims = 1:10)

# resolution range of 0.4-1.2 is usually good enough for 3k cells, larger dataset
# may need higher value for optimal resolution. increasing the resolution parameter
# value increases the number of clusters
ht67_cd205neg <- FindClusters(ht67_cd205neg, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(ht67_cd205neg), 5)


# ===========================================================


# run non-linear dimensional reduction
ht67_cd205neg <- RunUMAP(ht67_cd205neg, dims = 1:10)
# visualize (note that you can set `label = TRUE` or use the LabelClusters 
# function to help label individual clusters)
DimPlot(ht67_cd205neg, reduction = "umap", label = TRUE, label.size = 5)

# save output to RDS object
saveRDS(ht67_cd205neg, file = "ht67_cd205neg_tutorial_after_runumap.rds")
readRDS(ht67_cd205neg, file = "ht67_cd205neg_tutorial_after_runumap.rds")


# ===========================================================


# find differentially expressed genes (cluster biomarkers)
# find all cluster-defining markers (positive or negative) of a single cluster
# by default, it compares the cluster defined in ident.1 to all other cells
cluster2.markers <- FindMarkers(ht67_cd205neg, ident.1 = 2)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster_5_vs_0_3.markers <- FindMarkers(ht67_cd205neg, ident.1 = 5, ident.2 = c(0, 3))
head(cluster_5_vs_0_3.markers, n = 5)

# FindAllMarkers() automates this process for all clusters, but you can also test
# groups of clusters vs. each other, or against all cells.
FindAllMarkers(ht67_cd205neg)
# find all markers, return only positives, and save to object ht67_cd205neg.markers
ht67_cd205neg.markers <- FindAllMarkers(ht67_cd205neg, only.pos = TRUE)
# group the differential positive markers by cluster then subset those with
# fold change greater than 2 (i.e. avg_log2FC > 1)
ht67_cd205neg.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)


# ===========================================================


# differential gene expression visualization
# frequently used visualizations: VlnPlot(), FeaturePlot()
# other visualizations: RidgePlot(), CellScatter(), DotPlot()
VlnPlot(ht67_cd205neg, features = c("MS4A1", "CD79A"))
# you can plot raw counts as well
VlnPlot(ht67_cd205neg, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
# visualize gene expression on top of the umap plot
FeaturePlot(ht67_cd205neg, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))

# DoHeatmap() generates an expression heatmap for given cells and features. In
# this case, we are plotting the top 10 markers (or all markers if less than 10)
# for each cluster.
ht67_cd205neg.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(ht67_cd205neg, features = top10$gene) + NoLegend()


# ===========================================================

# save output to RDS object
saveRDS(ht67_cd205neg, file = "ht67_cd205neg_tutorial_before_assign_cell_type.rds")
readRDS(ht67_cd205neg, file = "ht67_cd205neg_tutorial_before_assign_cell_type.rds")


# assigning cell type identity to clusters
# add cell type names to an object, in order of the cluster ID from 0 to 8
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B Cell", "CD8 T",
                     "FCGR3A+ Mono", "NK", "DC", "Platelet")
# set the names of the cell types in new.cluster.ids as the corresponding levels
# of the clusters
names(new.cluster.ids) <- levels(ht67_cd205neg)
# assign cell type names to the corresponding clusters 0-8
ht67_cd205neg <- RenameIdents(ht67_cd205neg, new.cluster.ids)
# plot umap with cell type labels
DimPlot(ht67_cd205neg, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

library(ggplot2)
plot <- DimPlot(ht67_cd205neg, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
ggsave(filename = "output/images/ht67_cd205neg3k_umap.jpg", height = 7, width = 12, plot = plot, quality = 50)

# library(ggplot2)
# plot <- DimPlot(ht67_cd205neg, reduction = "umap", label = TRUE, label.size = 4.5)
#   + xlab("UMAP 1") 
#   + ylab("UMAP 2") 
#   + theme(axis.title = element_text(size = 18), 
#           legend.text = element_text(size = 18)) 
#   + guides(colour = guide_legend(override.aes = list(size = 10)))
# ggsave(filename = "output/images/ht67_cd205neg3k_umap.jpg", height = 7, width = 12, plot = plot, quality = 50)

saveRDS(ht67_cd205neg, file = "ht67_cd205neg3k_final.rds")

sessionInfo()

