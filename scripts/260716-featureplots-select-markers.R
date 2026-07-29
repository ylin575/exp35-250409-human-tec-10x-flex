# =============================================================================
# Composite UMAP FeaturePlot for the Fig 4B revision gene list
# Uses scCustomize::FeaturePlot_scCustom with viridis_light_high palette
# Rasterized point layer + 600 dpi for print-quality output at small file size.
# =============================================================================

library(Seurat)
library(scCustomize)
library(ggplot2)
library(patchwork)
library(ggrastr)

tec <- readRDS("data/rds-objects/260501-after-annotation-1.rds")

# -----------------------------------------------------------------------------
# 1. Gene list (matches the dotplot for internal consistency)
# -----------------------------------------------------------------------------
genes <- c("PTPRC","PSMB11","PRSS16","LY75","COL4A5","COL4A6","ITGB4","CD200",
           "CLU","AIRE","FEZF2","CD24","ITGA6","BCAM","CD74","PECAM1")

missing <- setdiff(genes, rownames(tec))
if (length(missing) > 0) {
  warning("Not in object, dropping: ", paste(missing, collapse = ", "))
  genes <- setdiff(genes, missing)
}

# -----------------------------------------------------------------------------
# 2. Build individual FeaturePlot_scCustom panels, rasterize points per panel
# -----------------------------------------------------------------------------
make_panel <- function(gene) {
  p <- FeaturePlot_scCustom(
    seurat_object = tec,
    reduction     = "umap",
    features      = gene,
    colors_use    = viridis_light_high,
    order         = FALSE,
    pt.size       = 0.1,
    raster        = FALSE      # we rasterize via ggrastr for finer control
  )
  
  # Rasterize just the point geom; text and legend stay vector
  rasterise(p, layers = "Point", dpi = 600)
}

panels <- lapply(genes, make_panel)

# -----------------------------------------------------------------------------
# 3. Compose into one figure
# -----------------------------------------------------------------------------
NCOL <- 4
composite <- wrap_plots(panels, ncol = NCOL)

# -----------------------------------------------------------------------------
# 4. Save (dimensions scale to panel count)
# -----------------------------------------------------------------------------
n_rows <- ceiling(length(panels) / NCOL)
w <- NCOL   * 5.4      # inches per panel column (slightly wider for viridis bar)
h <- n_rows * 4.8      # inches per panel row

ggsave("Fig4B_featureplots.pdf", composite,
       width = w, height = h, units = "in")
ggsave("Fig4B_featureplots.png", composite,
       width = w, height = h, units = "in", dpi = 600)

# -----------------------------------------------------------------------------
# 5. OPTIONAL -- PTPRC split by sample.id (strongest defense of the
#    EpCAM+autofluorescence sorting strategy against Reviewer 2's CD205 concern)
# -----------------------------------------------------------------------------
# ptprc_split <- FeaturePlot_scCustom(
#   seurat_object = tec,
#   reduction     = "umap",
#   features      = "PTPRC",
#   colors_use    = viridis_light_high,
#   order         = TRUE,
#   pt.size       = 0.4,
#   split.by      = "sample.id",
#   raster        = FALSE
# ) |> rasterise(layers = "Point", dpi = 600)
#
# ggsave("PTPRC_by_sample.pdf", ptprc_split, width = 14, height = 3.8)