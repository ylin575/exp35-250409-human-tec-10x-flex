# exp35-human-tec-10x-flex

# Load Seurat
library(dplyr)
library(Seurat)
library(scCustomize)
library(ggplot2)
library(harmony)

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

# run harmony
samples <- RunHarmony(samples, group.by.vars = "donor", reduction = "pca",
                      dims.use = 1:18, reduction.save = "harmony")

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
samples <- FindNeighbors(samples, reduction = "harmony", dims = 1:18)
samples <- FindClusters(samples, resolution = 0.5)
samples <- RunUMAP(samples, reduction = "harmony", dims = 1:18)


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
FeaturePlot_scCustom(samples, reduction = "umap", features = "PSMB11", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "LY75", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
# mTEC lo
FeaturePlot_scCustom(samples, reduction = "umap", features = "CCL19", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "KRT15", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "KRT19", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "FN1", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "CLU", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
# mTEC hi
FeaturePlot_scCustom(samples, reduction = "umap", features = "AIRE", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "FEZF2", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
# Basement Membrane Collagen-enriched TEC or "BMC TEC"
FeaturePlot_scCustom(samples, reduction = "umap", features = "ITGB4", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "COL4A5", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "COL4A6", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "COL8A1", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
# Nurse
FeaturePlot_scCustom(samples, reduction = "umap", features = "TBATA", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "PRSS16", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "CD3E", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "PTPRC", colors_use = viridis_light_high, label = FALSE, label.size = label_size, order = order)
# Mimetic TEC
FeaturePlot_scCustom(samples, reduction = "umap", features = "CD24", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
# All TECs
FeaturePlot_scCustom(samples, reduction = "umap", features = "TP63", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = order)
# Endothelial
FeaturePlot_scCustom(samples, reduction = "umap", features = "PECAM1", colors_use = viridis_light_high, label = TRUE, label.size = label_size, order = F)
# Other markers
FeaturePlot_scCustom(samples, reduction = "umap", features = "LY75", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "ITGB4", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "CD200", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "ITGA6", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "BCAM", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "EPCAM", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "CD24", colors_use = viridis_light_high, order = order)
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
FeaturePlot_scCustom(samples, reduction = "umap", features = "KRT1", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "KRT2", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "KRT3", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "KRT4", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "KRT5", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "KRT6", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "KRT7", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "KRT8", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "KRT9", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "KRT10", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "KRT11", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "KRT12", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "KRT13", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "KRT14", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "KRT15", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "KRT16", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "KRT17", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "KRT18", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "KRT19", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "KRT20", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "KRT21", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "KRT22", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "KRT23", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "KRT24", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "KRT25", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "KRT26", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "KRT27", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "KRT28", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "KRT29", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "KRT30", colors_use = viridis_light_high, order = order)

FeaturePlot_scCustom(samples, reduction = "umap", features = "CLU", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "AIRE", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "FEZF2", colors_use = viridis_light_high, order = order)

# bautista immature TEC-1
FeaturePlot_scCustom(samples, reduction = "umap", features = "ACTB", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "JUNB", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "FOS", colors_use = viridis_light_high, order = order)

# bautista immature TEC-2
FeaturePlot_scCustom(samples, reduction = "umap", features = "IGFBP5", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "NNMT", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "MAOA", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "DPYS", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "FKBP5", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "GLUL", colors_use = viridis_light_high, order = order)

# bautista mtec lo
FeaturePlot_scCustom(samples, reduction = "umap", features = "GABRA5", colors_use = viridis_light_high, order = order)
FeaturePlot_scCustom(samples, reduction = "umap", features = "LYPD1", colors_use = viridis_light_high, order = order)

# save object before annotation
saveRDS(samples, file = "data/rds-objects/260501-before-annotation.rds")
samples <- readRDS(file = "data/rds-objects/260501-before-annotation.rds")


# annotate cell types
test <- samples
new.cluster.ids <- c("cTEC",
                     "cTEC",
                     "mTEC lo",
                     "Nurse",
                     "BMC TEC",
                     "BMC TEC",
                     "mTEC lo",
                     "mTEC hi",
                     "mTEC lo",
                     "mTEC lo",
                     "BMC+ TEC",
                     "Mimetic",
                     "Mimetic",
                     "Mimetic",
                     "Nurse",
                     "Mimetic",
                     "Nurse",
                     "mTEC lo",
                     "mTEC hi",
                     "Mimetic",
                     "Mimetic",
                     "Endothelial",
                     "cTEC",
                     "cTEC",
                     "mTEC lo",
                     "cTEC",
                     "mTEC lo")
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
DimPlot(test, reduction = "umap", label = FALSE, split.by = "sample.id",
        cols = c("#D8A767","#F47D2B","#D24B27","#3D3D3D","#7E1416","#0570B0",
                 "#89C75F")
)


# save object after annotation
saveRDS(test, file = "data/rds-objects/260501-after-annotation-1.rds")
tec <- readRDS("data/rds-objects/260501-after-annotation-1.rds")


# 260706 plot for revision JI
genes <- c("PTPRC","PSMB11","PRSS16","LY75","COL4A5","COL4A6","ITGB4","CD200",
           "CLU","AIRE","FEZF2","CD24","PECAM1")

tec2 <- tec
my_order <- c("Nurse","cTEC","BMC+ TEC","mTEC I","mTEC II","Mimetic","Endothelial")
Idents(tec2) <- factor(Idents(tec2), levels = my_order)
VlnPlot(tec2, ncol = 3, features = genes)



