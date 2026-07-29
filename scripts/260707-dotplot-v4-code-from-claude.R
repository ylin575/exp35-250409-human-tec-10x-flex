# =============================================================================
# Fig 4B revision: split dotplot of TEC subset markers
# Two gene blocks: (1) general identity/cTEC/mTEC markers, (2) BMC-TEC / ECM markers
# Seurat v5.3.0
# =============================================================================

library(Seurat)
library(ggplot2)
library(patchwork)

# -----------------------------------------------------------------------------
# 1. Load object and set identities
# -----------------------------------------------------------------------------
# Replace with your actual object path / variable
tec <- readRDS("data/rds-objects/260501-after-annotation-1.rds")

# Confirm cluster identities are set to your TEC subset labels, e.g.:
# Idents(tec) <- tec$tec_subset

# Order clusters left-to-right as you want them displayed.
# Adjust these labels to match your actual annotation strings exactly.
cluster_order <- c("Endothelial","Mimetic","mTEC II","mTEC I","BMC+ TEC","cTEC","Nurse")
# If you also want nurse/mimetic TEC included on this panel, add them here, e.g.:
# cluster_order <- c("BMC-TEC", "cTEC", "nurse", "mTECI", "mTECII", "mimetic")

Idents(tec) <- factor(Idents(tec), levels = cluster_order)

# -----------------------------------------------------------------------------
# 2. Define gene blocks (named list -> Seurat facets each block into its own panel)
# -----------------------------------------------------------------------------
# Block 1: general / cTEC / mTEC identity markers (mirrors top block in Hannah's figure)
# Block 2: BMC-TEC / basement-membrane / ECM markers (mirrors bottom block)
#
# NOTE: Hannah's combined list did not include BCAM, which Reviewer 2 explicitly
# paired with CD49f (=ITGA6) for a joint module-score analysis. BCAM is included
# here as a commented-out optional addition -- add it if you want the dotplot
# itself to visualize BCAM alongside ITGA6 (this does NOT replace the separate
# module-score analysis Reviewer 2 asked for -- see fig4b_supp_modulescore.R note below).

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

# Sanity check: flag any requested genes not present in the object before plotting
all_genes <- unlist(gene_blocks)
missing_genes <- setdiff(all_genes, rownames(tec))
if (length(missing_genes) > 0) {
  warning("The following genes are not found in the Seurat object and will be dropped: ",
          paste(missing_genes, collapse = ", "))
  gene_blocks <- lapply(gene_blocks, function(g) intersect(g, rownames(tec)))
}

# -----------------------------------------------------------------------------
# 3. Build the dotplot
# -----------------------------------------------------------------------------
# DotPlot accepts a named list directly and will facet by list name (replicates
# the old SplitDotPlotGG two-panel layout seen in the mockup).

p <- DotPlot(
  object     = tec,
  features   = gene_blocks,
  assay      = "RNA",      # adjust if your normalized data lives in a different assay
  dot.scale  = 6,
  scale      = TRUE,       # z-scores average expression across groups (standard)
  cols       = c("#F7F4B4", "#1B3B6F")  # light yellow -> dark navy, approximates mockup gradient
) +
  RotatedAxis() +
  facet_grid(~ factor(feature.groups, levels = names(gene_blocks)),
             scales = "free_x", space = "free_x") +
  theme(
    axis.title       = element_blank(),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text       = element_text(face = "bold", size = 9),
    panel.spacing    = unit(0.6, "lines"),
    axis.text.x      = element_text(size = 9),
    axis.text.y      = element_text(size = 9)
  )

print(p)

# -----------------------------------------------------------------------------
# 4. Save
# -----------------------------------------------------------------------------
ggsave(
  filename = "Fig4B_dotplot.pdf",
  plot     = p,
  width    = 9,
  height   = 5,
  units    = "in"
)

# ggsave(
#   filename = "Fig4B_dotplot.png",
#   plot     = p,
#   width    = 9,
#   height   = 5,
#   units    = "in",
#   dpi      = 300
# )

# =============================================================================
# OPTIONAL: if you want a tile/heatmap style exactly matching Hannah's mockup
# (color-only, no dot-size encoding), replace step 3 with the block below.
# =============================================================================

library(dplyr)
library(tidyr)

make_heatmap_df <- function(tec, features) {
  avg <- AverageExpression(tec, features = features, assay = "RNA")$RNA
  avg_df <- as.data.frame(avg) %>%
    tibble::rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "cluster", values_to = "avg_exp")
  avg_df
}

df_all <- make_heatmap_df(tec, gene_blocks[[1]]) %>% mutate(block = names(gene_blocks)[1])
df_all$gene <- factor(df_all$gene, levels = rev(unlist(gene_blocks)))
df_all$cluster <- factor(df_all$cluster, levels = cluster_order)

ggplot(df_all, aes(x = cluster, y = gene, fill = avg_exp)) +
  geom_tile(color = "white") +
  facet_grid(block ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradientn(colors = c("#F7F4B4", "#3E9A9C", "#0B1F4B"), name = "Avg. Exp.") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
