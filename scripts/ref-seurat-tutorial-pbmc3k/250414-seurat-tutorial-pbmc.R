#> standard workflow outline
#> 1. load packages
#> 2. read in the count matrix data
#> 3. create a seurat object from the count matrix
#> 4. select and filter cells based on QC metrics
#> 5. data normalization
#> 6. detection of highly variable features/gene-expression across cells
#> 7. scaling the data
#> 8. perform linear dimensional reduction
#> 9. determine the 'dimensionality' of the dataset
#> 10. cluster the cells
#> 11. run non-linear dimensional reduction
#> 12. finding differentially expressed features (cluster biomarkers)
#> 13. assigning cell type identity to clusters


# Load Packages
library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC data set
pbmc.data <- Read10X(data.dir = "data/filtered_gene_bc_matrices/hg19/")

# Initialize the Seurat object with the raw (non-normalized data)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3,
                           min.features = 200)

# return some basic properties of the object: number of features, samples, type
# of assay
pbmc

# return more detailed structural information on the dataset
str(pbmc)
#> some notes on the structure of the dataset
#> a seurat dataset contains many types of information, divided into so-called
#> 'slots'
#> slots can have sub-slots
#> access to each slot or sub-slot is accomplished through using the @ and the $
#> symbol, these symbols are used alternately for each level of sub-slots. for
#> example, to access slot 'z' deep in a dataset, do a@b$c@d$e@z

# lets examine a few gens in the first thirty cells
head(pbmc.data)
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]


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
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# QC metrics are stored in @meta.data slot
# Show QC metrics for the first 5 cells
head(pbmc@meta.data, 5)

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# select only cells that meet the following criteria
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


# ===========================================================


# data normalization
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)


# ===========================================================


# find highly variable genes across cells and use these for downstream procedures
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


# ===========================================================


# scale data needed prior to running PCA
# Shifts the expression of each gene, so that the mean expression across cells is 0
# Scales the expression of each gene, so that the variance across cells is 1
# This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
# The results of this are stored in pbmc[["RNA"]]$scale.data
# By default, only variable features are scaled.
# You can specify the features argument to scale additional features
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)


# ===========================================================


# perform linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Seurat provides several useful ways of visualizing both cells and features that 
# define the PCA, including VizDimReduction(), DimPlot(), and DimHeatmap()

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca") + NoLegend()
DimHeatmap(pbmc, dims = 1:5, cells = 6, balanced = TRUE)
ElbowPlot(pbmc)


# ===========================================================


# cluster cells
# dims is the number of PCA dimensions to use, based on examining pca results above
pbmc <- FindNeighbors(pbmc, dims = 1:10)

# resolution range of 0.4-1.2 is usually good enough for 3k cells, larger dataset
# may need higher value for optimal resolution. increasing the resolution parameter
# value increases the number of clusters
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)


# ===========================================================


# run non-linear dimensional reduction
pbmc <- RunUMAP(pbmc, dims = 1:10)
# visualize (note that you can set `label = TRUE` or use the LabelClusters 
# function to help label individual clusters)
DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 5)

# save output to RDS object
saveRDS(pbmc, file = "pbmc_tutorial_after_runumap.rds")
readRDS(pbmc, file = "pbmc_tutorial_after_runumap.rds")


# ===========================================================


# find differentially expressed genes (cluster biomarkers)
# find all cluster-defining markers (positive or negative) of a single cluster
# by default, it compares the cluster defined in ident.1 to all other cells
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster_5_vs_0_3.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3))
head(cluster_5_vs_0_3.markers, n = 5)

# FindAllMarkers() automates this process for all clusters, but you can also test
# groups of clusters vs. each other, or against all cells.
FindAllMarkers(pbmc)
# find all markers, return only positives, and save to object pbmc.markers
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
# group the differential positive markers by cluster then subset those with
# fold change greater than 2 (i.e. avg_log2FC > 1)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)


# ===========================================================


# differential gene expression visualization
# frequently used visualizations: VlnPlot(), FeaturePlot()
# other visualizations: RidgePlot(), CellScatter(), DotPlot()
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
# visualize gene expression on top of the umap plot
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))

# DoHeatmap() generates an expression heatmap for given cells and features. In
# this case, we are plotting the top 10 markers (or all markers if less than 10)
# for each cluster.
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()


# ===========================================================

# save output to RDS object
saveRDS(pbmc, file = "pbmc_tutorial_before_assign_cell_type.rds")
readRDS(pbmc, file = "pbmc_tutorial_before_assign_cell_type.rds")


# assigning cell type identity to clusters
# add cell type names to an object, in order of the cluster ID from 0 to 8
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B Cell", "CD8 T",
                     "FCGR3A+ Mono", "NK", "DC", "Platelet")
# set the names of the cell types in new.cluster.ids as the corresponding levels
# of the clusters
names(new.cluster.ids) <- levels(pbmc)
# assign cell type names to the corresponding clusters 0-8
pbmc <- RenameIdents(pbmc, new.cluster.ids)
# plot umap with cell type labels
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

library(ggplot2)
plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
ggsave(filename = "output/images/pbmc3k_umap.jpg", height = 7, width = 12, plot = plot, quality = 50)

# library(ggplot2)
# plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4.5)
#   + xlab("UMAP 1") 
#   + ylab("UMAP 2") 
#   + theme(axis.title = element_text(size = 18), 
#           legend.text = element_text(size = 18)) 
#   + guides(colour = guide_legend(override.aes = list(size = 10)))
# ggsave(filename = "output/images/pbmc3k_umap.jpg", height = 7, width = 12, plot = plot, quality = 50)

saveRDS(pbmc, file = "pbmc3k_final.rds")

sessionInfo()
















