# 231030 parse analysis

# Loading libraries and setting location paths
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(cowplot)

rm(list = ls())

# data_path <- "/volume-general/analysis/data/pbmc/"
# fig_path <- "/volume-general/analysis/figures/pbmc/seurat/"


# Below we've included a few convenience functions for saving images and 
# reading/writing your Seurat object to disk. When you're working with larger 
# datasets, it's usually a good idea to save your progress after computationally 
# intensive steps so you can back track if you wish to do so.

# Convenience functions
SaveFigure <- function(plots, name, type = "png", width, height, res){
  if(type == "png") {
    png(paste0(fig_path, name, ".", type),
        width = width, height = height, units = "in", res = 200)
  } else {
    pdf(paste0(fig_path, name, ".", type),
        width = width, height = height)
  }
  print(plots)
  dev.off()
}

SaveObject <- function(object, name){
  saveRDS(object, paste0(data_path, name, ".RDS"))
}

ReadObject <- function(name){
  readRDS(paste0(data_path, name, ".RDS"))
}


# Reading in data

# After reading in the data we'll perform basic filtering a on our expression 
# matrix to remove low quality cells and uninformative genes. The parameter 
# "min_genes" will keep cells that have at least 300 genes, and similarly, 
# "min_cells" will keep genes that are expressed in at least 5 cells. 
# Note: Seurat version 4.1 includes a convenience function to read Parse data from 
# the DGE folder. If you would like to use this function, please skip the code
# block below and see the section "Reading in data with Seurat >= 4.1"

mat_path <- "rawdata/all-sample/DGE_filtered"
mat <- ReadParseBio(mat_path)

# Check to see if empty gene names are present, add name if name is absent
table(rownames(mat) == "")
rownames(mat)[rownames(mat) == ""] <- "unknown"

# Read in cell meta data
cell_meta <- read.csv(paste0(mat_path, "/cell_metadata.csv"), row.names = 1)

# Create object
hthy <- CreateSeuratObject(mat, min_genes = 100, min_cells = 100,
                           names.feild = 0, meta.data = cell_meta)


# When we create our Seurat object the plate well numbers (column names in the 
# expression matrix) from the experiment will automatically be assigned to the 
# cell identity slot. In other words, the program assumes this is how we want to 
# initially classify our cells. In general, we would like to avoid this behavior 
# so there isn't a bias towards a particular cell class when removing outliers.

# Setting our initial cell class to a single type, this will change after 
# clustering. 
hthy@meta.data$orig.ident <- factor(rep("hthy", nrow(hthy@meta.data)))
Idents(hthy) <- hthy@meta.data$orig.ident

# copied below 2 lines from parse tutorial, doesn't work (can't find functions)
SaveObject(hthy, "seurat_obj_before_QC")
hthy <- ReadObject("seurat_obj_before_QC")

# use saveRSD instead
saveRDS(hthy, file = "data/rds_objects/seurat_obj_before_QC_231120.rds")

#Cell QC

hthy[["percent.mt"]] <- PercentageFeatureSet(hthy, pattern = "^MT-")
plot <- VlnPlot(hthy, pt.size = 0.10,
                features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                ncol = 3)
plot
# SaveFigure(plot, "vln_QC", width = 12, height = 6)

# plot1 <- FeatureScatter(hthy, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(hthy, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# SaveFigure((plot1 + plot2),"scatter_QC", width = 12, height = 6, res = 200)

# feature scatter
plot1 <- FeatureScatter(hthy, feature1 = "nCount_RNA", feature2 = "percent.mt",
                        group.by = "sample")
plot2 <- FeatureScatter(hthy, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                        group.by = "sample")

plot1 + plot2
# Visualize QC metrics as histograms, find best cut off point

# nFeature_RNA
hist(hthy@meta.data$nFeature_RNA[hthy@meta.data$sample=="Ht9"], 
     breaks=500, xlim=c(0,5000),
     xlab = "nFeature_RNA", ylab = "Frequency", main = "HT9")
hist(hthy@meta.data$nFeature_RNA[hthy@meta.data$sample=="HT11"], 
     breaks=500, xlim=c(0,5000),
     xlab = "nFeature_RNA", ylab = "Frequency", main = "HT11")
hist(hthy@meta.data$nFeature_RNA[hthy@meta.data$sample=="Ht12"], 
     breaks=500, xlim=c(0,5000),
     xlab = "nFeature_RNA", ylab = "Frequency", main = "HT12")
hist(hthy@meta.data$nFeature_RNA[hthy@meta.data$sample=="Ht14"], 
     breaks=300, xlim=c(0,5000),
     xlab = "nFeature_RNA", ylab = "Frequency", main = "HT14")

# nCount_RNA
hist(hthy@meta.data$nCount_RNA[hthy@meta.data$sample=="Ht9"], 
     breaks=500, xlim=c(0,20000),
     xlab = "nCount_RNA", ylab = "Frequency", main = "HT9")
hist(hthy@meta.data$nCount_RNA[hthy@meta.data$sample=="HT11"], 
     breaks=500, xlim=c(0,20000),
     xlab = "nCount_RNA", ylab = "Frequency", main = "HT11")
hist(hthy@meta.data$nCount_RNA[hthy@meta.data$sample=="Ht12"], 
     breaks=500, xlim=c(0,20000),
     xlab = "nCount_RNA", ylab = "Frequency", main = "HT12")
hist(hthy@meta.data$nCount_RNA[hthy@meta.data$sample=="Ht14"], 
     breaks=500, xlim=c(0,20000),
     xlab = "nCount_RNA", ylab = "Frequency", main = "HT14")

# percent.mt
hist(hthy@meta.data$percent.mt[hthy@meta.data$sample=="Ht9"], 
     breaks=200, xlim=c(0,60),
     xlab = "percent.mt", ylab = "Frequency", main = "HT9")
hist(hthy@meta.data$percent.mt[hthy@meta.data$sample=="HT11"], 
     breaks=200, xlim=c(0,60),
     xlab = "percent.mt", ylab = "Frequency", main = "HT11")
hist(hthy@meta.data$percent.mt[hthy@meta.data$sample=="Ht12"], 
     breaks=200, xlim=c(0,60),
     xlab = "percent.mt", ylab = "Frequency", main = "HT12")
hist(hthy@meta.data$percent.mt[hthy@meta.data$sample=="Ht14"], 
     breaks=200, xlim=c(0,60),
     xlab = "percent.mt", ylab = "Frequency", main = "HT14")

# Compare QC metrics across samples using boxplot
boxplot(hthy@meta.data$nCount_RNA ~ hthy@meta.data$sample, 
          xlab = "sample #", ylab = "nCount_RNA")
boxplot(hthy@meta.data$nFeature_RNA ~ hthy@meta.data$sample, 
        xlab = "sample #", ylab = "nFeature_RNA")
boxplot(hthy@meta.data$percent.mt ~ hthy@meta.data$sample, 
        xlab = "sample #", ylab = "percent.mt")

#violin plot in linear scale
VlnPlot(hthy, pt.size = 0.001, group.by = "sample",
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)
#violin plot in log scale
VlnPlot(hthy, pt.size = 0.001, group.by = "sample", log = TRUE,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)

