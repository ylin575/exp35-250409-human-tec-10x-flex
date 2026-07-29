# =============================================================================
# Fig 4B revision: split dotplot of TEC subset markers
# Produces TWO versions:
#   (A) SCALED   -> color = z-score of log1p(avg.exp), per gene   [Seurat default]
#   (B) UNSCALED -> color = log1p(avg.exp), comparable across genes ("absolute")
#
# Verified against Seurat v5.3.0 source (satijalab/seurat tag v5.3.0):
#   avg.exp          = mean(expm1(x))              # linear-scale mean per cluster
#   scale = TRUE  -> scale(log1p(avg.exp)), clipped to [col.min, col.max]=[-2.5,2.5]
#   scale = FALSE -> log1p(avg.exp)                # non-negative (NOT natural log)
#   dot size = percent expressing, in BOTH cases
#   col.min/col.max are inert when scale = FALSE
# =============================================================================

library(Seurat)
library(ggplot2)
library(patchwork)

# -----------------------------------------------------------------------------
# 1. Load object and set identities
# -----------------------------------------------------------------------------
tec <- readRDS("data/rds-objects/260501-after-annotation-1.rds")
# Idents(tec) <- tec$tec_subset   # ensure identities are your TEC subset labels

cluster_order <- c("Endothelial","Mimetic","mTEC II","mTEC I","BMC+ TEC","cTEC","Nurse")
# add nurse/mimetic if you want them on this panel, e.g.:
# cluster_order <- c("BMC-TEC", "cTEC", "nurse", "mTECI", "mTECII", "mimetic")

Idents(tec) <- factor(Idents(tec), levels = cluster_order)

# -----------------------------------------------------------------------------
# 2. Gene blocks (named list -> Seurat facets each block into its own panel)
# -----------------------------------------------------------------------------
# NOTE: BCAM (paired with ITGA6=CD49f by Reviewer 2 for the module score) is
# commented out. Adding it here visualizes it on the dotplot but does NOT
# substitute for the separate AddModuleScore() analysis Reviewer 2 requested.

gene_blocks <- list(
  "Select Genes" = c(
    "PTPRC",   # Ptprc -- reviewer 2 explicit request
    "LY75",
    "PSMB11",
    "PRSS16",
    "COL4A5",
    "COL4A6",
    "COL8A1",
    "ITGB4",
    "CD200",
    "CCL19",
    "CLU",
    "KRT15",
    "AIRE",
    "FEZF2",
    "CD24",
    "ITGA6",   # = CD49f
    "BCAM",  # <- uncomment to add; paired with ITGA6 per Reviewer 2's module request
    "CD74",
    "PECAM1"
  )
)

# Drop any genes absent from the object (with a warning) before plotting
all_genes <- unlist(gene_blocks)
missing_genes <- setdiff(all_genes, rownames(tec))
if (length(missing_genes) > 0) {
  warning("Genes not found and dropped: ", paste(missing_genes, collapse = ", "))
  gene_blocks <- lapply(gene_blocks, function(g) intersect(g, rownames(tec)))
}

# -----------------------------------------------------------------------------
# 3. Helper to build a styled dotplot (scaled or unscaled)
# -----------------------------------------------------------------------------
make_dotplot <- function(obj, features, scale_expr, cols_pair, title_txt) {
  p <- DotPlot(
    object    = obj,
    features  = features,
    assay     = "RNA",         # adjust if your normalized data is elsewhere
    dot.scale = 6,
    scale     = scale_expr,
    cols      = cols_pair
  ) +
    RotatedAxis() +
    facet_grid(~ factor(feature.groups, levels = names(features)),
               scales = "free_x", space = "free_x") +
    labs(title = title_txt) +
    theme(
      axis.title       = element_blank(),
      plot.title       = element_text(size = 10, face = "bold"),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text       = element_text(face = "bold", size = 9),
      panel.spacing    = unit(0.6, "lines"),
      axis.text.x      = element_text(size = 9),
      axis.text.y      = element_text(size = 9)
    )
  return(p)
}

# -----------------------------------------------------------------------------
# 4. (A) SCALED version  (color = z-scored log1p avg expression, per gene)
# -----------------------------------------------------------------------------
p_scaled <- make_dotplot(
  obj        = tec,
  features   = gene_blocks,
  scale_expr = TRUE,
  cols_pair  = c("#F7F4B4", "#1B3B6F"),   # light yellow -> dark navy
  title_txt  = "Scaled expression (z-score of log1p avg, per gene)"
)

# -----------------------------------------------------------------------------
# 5. (B) UNSCALED / ABSOLUTE version  (color = log1p(avg.exp), cross-gene comparable)
# -----------------------------------------------------------------------------
# Recommended for defending "BMC cells express little/no CD200" (Reviewer 2),
# because it shows true magnitude rather than per-gene relative color.
p_absolute <- make_dotplot(
  obj        = tec,
  features   = gene_blocks,
  scale_expr = FALSE,
  cols_pair  = c("lightgrey", "#1B3B6F"), # grey (low) -> navy (high), non-negative range
  title_txt  = "scale_expr = FALSE"
)

print(p_scaled)
print(p_absolute)

# -----------------------------------------------------------------------------
# 6. Save both
# -----------------------------------------------------------------------------
ggsave("Fig4B_dotplot_scaled.pdf",   p_scaled,   width = 9, height = 5, units = "in")
ggsave("Fig4B_dotplot_absolute.pdf", p_absolute, width = 9, height = 5, units = "in")
# ggsave("Fig4B_dotplot_scaled.png",   p_scaled,   width = 9, height = 5, units = "in", dpi = 300)
# ggsave("Fig4B_dotplot_absolute.png", p_absolute, width = 9, height = 5, units = "in", dpi = 300)

# Optional side-by-side for quick comparison:
# ggsave("Fig4B_dotplot_compare.pdf", p_scaled / p_absolute, width = 9, height = 10, units = "in")

# =============================================================================
# NOTE on matching Hannah's mockup exactly:
# The mockup shows a SEPARATE color scale per gene block (e.g. ~0.5-2 for the
# identity block, ~0.2-0.6 for the ECM block). DotPlot with a named list uses
# ONE shared color scale across both facets. To reproduce independent per-block
# scales, build two separate DotPlots (one per block) and combine with
# patchwork (p_block1 | p_block2). Ask if you want that variant.
# =============================================================================