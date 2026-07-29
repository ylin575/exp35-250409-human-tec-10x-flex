# =============================================================================
# Composite UMAP FeaturePlot for the Fig 4B revision gene list
# Directly answers Reviewer 2's request: "The authors don't show a Ptprc
# feature/violin plot in their scRNAseq dataset. They should add it..."
#
# Uses the same gene list as the dotplot for internal consistency, and
# rasterizes the point layer to keep the PDF small (large UMAPs of thousands
# of cells otherwise blow up to tens of MB per figure).
# =============================================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(ggrastr)

tec <- readRDS("data/rds-objects/260501-after-annotation-1.rds")

# -----------------------------------------------------------------------------
# 1. Gene list (matches the dotplot; adjust BCAM in/out per your final decision)
# -----------------------------------------------------------------------------
genes <- c("PTPRC","PSMB11","PRSS16","LY75","COL4A5","COL4A6","ITGB4","CD200",
             "CLU","AIRE","FEZF2","CD24","ITGA6","BCAM","CD74","PECAM1")

# Drop genes not present in the object
missing <- setdiff(genes, rownames(tec))
if (length(missing) > 0) {
  warning("Not in object, dropping: ", paste(missing, collapse = ", "))
  genes <- setdiff(genes, missing)
}

# -----------------------------------------------------------------------------
# 2. Build individual FeaturePlots, rasterize the point layer per panel
# -----------------------------------------------------------------------------
# Building panels one at a time (rather than one FeaturePlot() call with all
# genes) gives us clean control over rasterization and consistent theming.

make_panel <- function(gene) {
  p <- FeaturePlot(
    tec,
    features    = gene,
    order       = TRUE,        # plot high-expression cells on top
    pt.size     = 0.4,
    raster      = FALSE,       # we rasterize ourselves via ggrastr for control
    combine     = TRUE
  ) +
    scale_color_gradientn(colors = c("lightgrey", "#1B3B6F")) +
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

panels <- lapply(genes, make_panel)

# -----------------------------------------------------------------------------
# 3. Compose into one figure
# -----------------------------------------------------------------------------
# Layout: 4 columns is a good default for 15-18 genes -> 4-5 rows
NCOL <- 4
composite <- wrap_plots(panels, ncol = NCOL)

# -----------------------------------------------------------------------------
# 4. Save (dimensions scale to panel count)
# -----------------------------------------------------------------------------
n_rows <- ceiling(length(panels) / NCOL)
w <- NCOL   * 2.5      # inches per panel column
h <- n_rows * 2.2      # inches per panel row

ggsave("Fig4B_featureplots.pdf", composite,
       width = w, height = h, units = "in")
ggsave("Fig4B_featureplots.png", composite,
       width = w, height = h, units = "in", dpi = 300)

# -----------------------------------------------------------------------------
# 5. OPTIONAL -- side-by-side per sort condition (paired with your existing
#    UMAP split-by-sort figure)
# -----------------------------------------------------------------------------
# If you want to show that, e.g., PTPRC is specifically absent in sorted CD205+
# TECs (which directly addresses Reviewer 2's concern about the EpCAM+auto-
# fluorescence sorting strategy), split-by sample.id:
#
# make_split_panel <- function(gene) {
#   p <- FeaturePlot(tec, features = gene, split.by = "sample.id",
#                    order = TRUE, pt.size = 0.4, raster = FALSE) &
#     scale_color_gradientn(colors = c("lightgrey", "#1B3B6F")) &
#     theme_void(base_size = 9) &
#     theme(plot.title = element_text(face = "italic", size = 10, hjust = 0.5))
#   rasterise(p, layers = "Point", dpi = 300)
# }
#
# ptprc_split <- make_split_panel("PTPRC")
# ggsave("PTPRC_by_sample.pdf", ptprc_split, width = 12, height = 3.5)