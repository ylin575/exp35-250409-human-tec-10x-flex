# =============================================================================
# Fig 4B (rectangular / heatmap style, matching Hannah's mockup layout)
# RELATIVE scale version: per-gene (row) z-score across clusters
#
# Contrast with the absolute version:
#   AverageExpression default -> mean(expm1(data)) = ABSOLUTE linear expression
#   This script                -> scale(log1p(avg)) per gene = RELATIVE (z-score)
#
# The log1p-then-scale transform matches DotPlot(scale = TRUE) exactly, so this
# heatmap is internally consistent with the scaled dotplot you already generated.
# Verified against Seurat v5.3.0 source (satijalab/seurat tag v5.3.0).
# =============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

# -----------------------------------------------------------------------------
# 1. Object + identities (same as the dotplot script)
# -----------------------------------------------------------------------------
tec <- readRDS("data/rds-objects/260501-after-annotation-1.rds")
# Idents(tec) <- tec$tec_subset

cluster_order <- c("Nurse", "cTEC", "BMC+ TEC", "mTEC I", "mTEC II", "Mimetic", "Endothelial")
# cluster_order <- c("Endothelial","Mimetic","mTEC II","mTEC I","BMC+ TEC","cTEC","Nurse")
Idents(tec) <- factor(Idents(tec), levels = cluster_order)

# -----------------------------------------------------------------------------
# 2. Gene blocks (same two-block structure as the mockup)
# -----------------------------------------------------------------------------
gene_blocks <- list(
  "Select Genes" = c(
    "PTPRC",   # Ptprc -- reviewer 2 explicit request
    "LY75",
    "PSMB11",
    "PRSS16",
    "CCL25",
    "COL4A5",
    "COL4A6",
    "COL8A1",
    "IGFBP5",
    "CDH13",
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
    "PECAM1",
    "EPCAM",
    "TP63",
    "DLK2",
    "KRT14",
    "FOXN1",
    "ACTA2",
    "CLDN4"
  )
)

# bautista immature tec markers
gene_blocks <- list(
  "Select Genes" = c(
    "IGFBP6",
    "KRT15",
    "ZBED2",
    "CDH13",
    "CCN2", # note that CCN2 is the new name for CTGF used in Bautista paper
    "IGFBP5",
    "NNMT",
    "MAOA",
    "DPYS",
    "FKBP5",
    "GLUL"
  )
)


# bautista mtec lo markers
gene_blocks <- list(
  "Select Genes" = c(
    "GABRA5",
    "LYPD1"
  )
)


# bautista mtec hi markers
gene_blocks <- list(
  "Select Genes" = c(
    "CLEC7A",
    "MARCO",
    "FXYD2",
    "FXYD3",
    "IL4I1",
    "CHI3L1",
    "CD70",
    "TNFRSF9"
  )
)

# bautista mtec hi markers
gene_blocks <- list(
  "Select Genes" = c(
    "CLEC7A",
    "MARCO",
    "FXYD2",
    "FXYD3",
    "IL4I1",
    "CHI3L1",
    "CD70",
    "TNFRSF9"
  )
)

# bautista corneocyte-like mtec markers
gene_blocks <- list(
  "Select Genes" = c(
    "FXYD3",
    "IL1RN",
    "LYPD2"
  )
)

# bautista genes in immature tecs that increase expression with age
gene_blocks <- list(
  "Select Genes" = c(
    "IGFBP5",
    "IGKC",
    "NNMT",
    "TXNIP",
    "SAA1",
    "MAOA",
    "FAM107A",
    "TSC22D3",
    "CSTA",
    "FKBP5",
    "DPYS",
    "ZBED2"
  )
)


all_genes <- unlist(gene_blocks)
missing_genes <- setdiff(all_genes, rownames(tec))
if (length(missing_genes) > 0) {
  warning("Genes not found and dropped: ", paste(missing_genes, collapse = ", "))
  gene_blocks <- lapply(gene_blocks, function(g) intersect(g, rownames(tec)))
  all_genes <- unlist(gene_blocks)
}

# -----------------------------------------------------------------------------
# 3. Compute per-gene relative (z-scored) expression
# -----------------------------------------------------------------------------
# (a) absolute linear average from Seurat (genes x clusters)
avg_mat <- AverageExpression(tec, features = all_genes, assay = "RNA")$RNA
avg_mat <- avg_mat[all_genes, cluster_order, drop = FALSE]   # enforce order

# (b) log1p to match the DotPlot(scale=TRUE) convention, THEN row z-score.
#     scale() works column-wise, so transpose -> scale -> transpose back to
#     z-score each GENE (row) across clusters.
log_mat <- log1p(avg_mat)
z_mat   <- t(scale(t(log_mat)))

# (c) guard: genes with zero variance across clusters produce NaN -> set to 0
z_mat[is.nan(z_mat)] <- 0

# If you'd rather z-score the RAW linear average (skip log1p), replace (b) with:
#   z_mat <- t(scale(t(avg_mat))); z_mat[is.nan(z_mat)] <- 0

# -----------------------------------------------------------------------------
# 4. Long format for ggplot, tagged by gene block
# -----------------------------------------------------------------------------
block_of <- stack(lapply(gene_blocks, identity))          # gene -> block map
block_lookup <- setNames(as.character(block_of$ind), block_of$values)

df <- as.data.frame(z_mat) %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "cluster", values_to = "z") %>%
  mutate(
    block   = factor(block_lookup[gene], levels = names(gene_blocks)),
    gene    = factor(gene, levels = rev(all_genes)),        # top-to-bottom order
    cluster = factor(cluster, levels = cluster_order)
  )

# -----------------------------------------------------------------------------
# 5. Plot (diverging fill centered at 0 -- appropriate for z-scores)
# -----------------------------------------------------------------------------
p <- ggplot(df, aes(x = cluster, y = gene, fill = z)) +
  geom_tile(color = "white", linewidth = 0.4) +
  facet_grid(block ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(
    low = "#F7E27A", mid = "grey95", high = "#1B3B6F",
    midpoint = 0, name = "Relative\nexpr. (z)"
  ) +
  labs(title = "Select Markers") +
  theme_minimal(base_size = 10) +
  theme(
    axis.title       = element_blank(),
    plot.title       = element_text(size = 10, face = "bold"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    strip.text.y     = element_blank(),
    panel.grid       = element_blank(),
    panel.spacing    = unit(0.4, "lines")
  )

print(p)

ggsave("Fig4B_heatmap_relative.pdf", p, width = 6, height = 6, units = "in")
ggsave("Fig4B_heatmap_relative.png", p, width = 6, height = 6, units = "in", dpi = 300)