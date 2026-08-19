# need to load the 10x flex scrnaseq dataset from 260807-fresh-start script

# -----------------------------------------------------------------------------
# nclTEC module score -- whole dataset, across all annotated cell types
# (self-contained -- re-derives everything from scratch)
# -----------------------------------------------------------------------------
nclTEC_genes <- c("CCL25", "PRSS16", "PSMB11", "NEURL2", "CBX4", "FOXN1",
                  "BMP8A", "TP63", "KRT5", "KRT8", "JUNB")

missing_genes <- setdiff(nclTEC_genes, rownames(samples))
if (length(missing_genes) > 0) {
  warning("Not in object, dropping: ", paste(missing_genes, collapse = ", "))
  nclTEC_genes <- setdiff(nclTEC_genes, missing_genes)
}
stopifnot(length(nclTEC_genes) > 0)

# sanity check: confirm cell_type_unintegrated levels before scoring
cat("Current cell_type_unintegrated levels:\n")
print(table(samples$cell_type_unintegrated))

samples <- AddModuleScore(
  samples, features = list(nclTEC_genes),
  name = "nclTEC_score_raw", verbose = TRUE
)
samples$nclTEC_score <- samples$nclTEC_score_raw1
samples$nclTEC_score_raw1 <- NULL

# sanity check: confirm score column created with no NAs
stopifnot("nclTEC_score" %in% colnames(samples@meta.data))
stopifnot(!any(is.na(samples$nclTEC_score)))

cat("\n=== nclTEC score by cell type (whole dataset) ===\n")
nclTEC_summary <- samples@meta.data %>%
  group_by(cell_type_unintegrated) %>%
  summarise(n = n(), mean_score = round(mean(nclTEC_score), 3),
            median_score = round(median(nclTEC_score), 3), .groups = "drop") %>%
  arrange(desc(mean_score))
print(nclTEC_summary)

celltype_order <- c("Nurse", "cTEC", "BMC TEC", "mTEC lo", "mTEC hi",
                    "Mimetic", "Bipolar TEC", "Endothelial", "Pericyte", "Erythroid")

# sanity check: confirm celltype_order matches the object's actual current levels
current_levels <- levels(droplevels(samples$cell_type_unintegrated))
stopifnot(setequal(celltype_order, current_levels))

samples$cell_type_unintegrated <- factor(samples$cell_type_unintegrated, levels = celltype_order)

p_nclTEC_whole <- VlnPlot(samples, features = "nclTEC_score",
                          group.by = "cell_type_unintegrated", pt.size = 0) +
  labs(title = "nclTEC module score across whole dataset", x = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_nclTEC_whole)
ggsave(paste0(DATE_PREFIX, "-wholedataset-nclTEC-score-violin.pdf"),
       p_nclTEC_whole, width = 8, height = 5)

write.csv(nclTEC_summary, paste0(DATE_PREFIX, "-wholedataset-nclTEC-score-summary.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# nclTEC module score FeaturePlot
# -----------------------------------------------------------------------------
p_umap <- FeaturePlot(samples, features = "nclTEC_score",
                      reduction = "umap_unintegrated", order = TRUE) +
  scale_color_viridis_c(option = "viridis") +
  labs(title = "nclTEC module score, UMAP (unintegrated)")
print(p_umap)
ggsave(paste0(DATE_PREFIX, "-PolyKRT-module-score-umap-unintegrated.pdf"),
       p_umap, width = 6, height = 5)

# -----------------------------------------------------------------------------
# 2D scatter: nclTEC score (x) vs. PolyKRT score (y), whole dataset
# (self-contained -- re-derives both scores from scratch)
# -----------------------------------------------------------------------------
nclTEC_genes <- c("CCL25", "PRSS16", "PSMB11", "NEURL2", "CBX4", "FOXN1",
                  "BMP8A", "TP63", "KRT5", "KRT8", "JUNB")
polyKRT_module_genes <- c(
  "CEBPD", "CLU", "FN1", "IFITM3", "TIMP1", "VCAM1", "TAGLN",
  "BCAM", "LIFR", "CH25H", "CCN2", "CCL19",
  "KRT13", "KRT14", "KRT15", "KRT17", "KRT19"
)

nclTEC_genes <- intersect(nclTEC_genes, rownames(samples))
polyKRT_module_genes <- intersect(polyKRT_module_genes, rownames(samples))
stopifnot(length(nclTEC_genes) > 0, length(polyKRT_module_genes) > 0)

samples <- AddModuleScore(samples, features = list(nclTEC_genes),
                          name = "nclTEC_scatter_raw", verbose = TRUE)
samples <- AddModuleScore(samples, features = list(polyKRT_module_genes),
                          name = "PolyKRT_scatter_raw", verbose = TRUE)

samples$nclTEC_score_2d  <- samples$nclTEC_scatter_raw1
samples$PolyKRT_score_2d <- samples$PolyKRT_scatter_raw1
samples$nclTEC_scatter_raw1 <- NULL
samples$PolyKRT_scatter_raw1 <- NULL

# sanity check: confirm both score columns exist with no NAs
stopifnot(all(c("nclTEC_score_2d", "PolyKRT_score_2d") %in% colnames(samples@meta.data)))
stopifnot(!any(is.na(samples$nclTEC_score_2d)), !any(is.na(samples$PolyKRT_score_2d)))

scatter_df <- FetchData(samples, vars = c("nclTEC_score_2d", "PolyKRT_score_2d",
                                          "cell_type_unintegrated"))

p_scatter_2d <- ggplot(scatter_df, aes(x = nclTEC_score_2d, y = PolyKRT_score_2d)) +
  geom_point(aes(color = cell_type_unintegrated), size = 0.6, alpha = 0.85) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  labs(title = "nclTEC score vs. PolyKRT score, whole dataset",
       x = "nclTEC module score", y = "PolyKRT module score", color = "Cell type") +
  theme_minimal(base_size = 11)

print(p_scatter_2d)
ggsave(paste0(DATE_PREFIX, "-nclTEC-vs-PolyKRT-scatter-wholedataset.pdf"),
       p_scatter_2d, width = 9, height = 7)


# -----------------------------------------------------------------------------
# Add BMC TEC highlight overlay on top of the existing scatter
# -----------------------------------------------------------------------------
p_scatter_2d <- ggplot(scatter_df, aes(x = nclTEC_score_2d, y = PolyKRT_score_2d)) +
  geom_point(aes(color = cell_type_unintegrated), size = 0.6, alpha = 0.85) +
  geom_point(data = filter(scatter_df, cell_type_unintegrated == "BMC TEC"),
             aes(x = nclTEC_score_2d, y = PolyKRT_score_2d),
             color = "black", size = 1.2, alpha = 0.9, shape = 21, fill = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  labs(title = "nclTEC score vs. PolyKRT score, whole dataset (BMC TEC highlighted)",
       x = "nclTEC module score", y = "PolyKRT module score", color = "Cell type") +
  theme_minimal(base_size = 11)

print(p_scatter_2d)
ggsave(paste0(DATE_PREFIX, "-nclTEC-vs-PolyKRT-scatter-wholedataset-BMCTEC-highlight.pdf"),
       p_scatter_2d, width = 9, height = 7)



# -----------------------------------------------------------------------------
# Save one scatter plot per cell type, each highlighting a different type,
# with an automatically generated matching title
# -----------------------------------------------------------------------------
celltype_order <- c("Nurse", "cTEC", "BMC TEC", "mTEC lo", "mTEC hi",
                    "Mimetic", "Bipolar TEC", "Endothelial", "Pericyte", "Erythroid")

# sanity check: confirm celltype_order matches the object's actual current levels
stopifnot(setequal(celltype_order, levels(droplevels(samples$cell_type_unintegrated))))

for (ct in celltype_order) {
  p_ct_highlight <- ggplot(scatter_df, aes(x = nclTEC_score_2d, y = PolyKRT_score_2d)) +
    geom_point(aes(color = cell_type_unintegrated), size = 0.6, alpha = 0.85) +
    geom_point(data = filter(scatter_df, cell_type_unintegrated == ct),
               aes(x = nclTEC_score_2d, y = PolyKRT_score_2d),
               color = "black", size = 1.2, alpha = 0.9, shape = 21, fill = "red") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
    labs(title = paste0("nclTEC score vs. PolyKRT score, whole dataset (", ct, " highlighted)"),
         x = "nclTEC module score", y = "PolyKRT module score", color = "Cell type") +
    theme_minimal(base_size = 11)
  
  ggsave(paste0(DATE_PREFIX, "-nclTEC-vs-PolyKRT-scatter-", gsub(" ", "", ct), "-highlight.pdf"),
         p_ct_highlight, width = 9, height = 7)
}


# -----------------------------------------------------------------------------
# Rotated oval (elliptical) gate -- adds a rotation angle to the point-in-
# ellipse test and the plotted outline
# -----------------------------------------------------------------------------
center_x <- 1.5
center_y <- 0.4
radius_x <- 0.8
radius_y <- 0.9
angle_deg <- 25   # counter-clockwise rotation, in degrees
angle_rad <- angle_deg * pi / 180

# Point-in-rotated-ellipse test: translate to center, then rotate the point
# by -angle (undoing the ellipse's rotation) before applying the standard
# axis-aligned ellipse formula
scatter_df_ellipse <- scatter_df %>%
  tibble::rownames_to_column("cell") %>%
  mutate(
    x_centered = nclTEC_score_2d - center_x,
    y_centered = PolyKRT_score_2d - center_y,
    x_rot = x_centered * cos(-angle_rad) - y_centered * sin(-angle_rad),
    y_rot = x_centered * sin(-angle_rad) + y_centered * cos(-angle_rad),
    ellipse_dist = (x_rot / radius_x)^2 + (y_rot / radius_y)^2,
    in_ellipse = ellipse_dist <= 1
  )

# sanity check: confirm a non-trivial, non-total subset was captured
n_in_ellipse <- sum(scatter_df_ellipse$in_ellipse)
cat("Cells within rotated ellipse gate:", n_in_ellipse, "of", nrow(scatter_df_ellipse), "\n")
stopifnot(n_in_ellipse > 0, n_in_ellipse < nrow(scatter_df_ellipse))

cat("\nCell type composition within the ellipse:\n")
print(table(samples$cell_type_unintegrated[scatter_df_ellipse$cell[scatter_df_ellipse$in_ellipse]]))

# -----------------------------------------------------------------------------
# Build rotated ellipse outline coordinates for plotting
# -----------------------------------------------------------------------------
theta <- seq(0, 2 * pi, length.out = 200)
ellipse_outline <- data.frame(
  x_unrot = radius_x * cos(theta),
  y_unrot = radius_y * sin(theta)
) %>%
  mutate(
    x = center_x + x_unrot * cos(angle_rad) - y_unrot * sin(angle_rad),
    y = center_y + x_unrot * sin(angle_rad) + y_unrot * cos(angle_rad)
  )

p_ellipse_check <- ggplot(scatter_df_ellipse, aes(x = nclTEC_score_2d, y = PolyKRT_score_2d)) +
  geom_point(aes(color = in_ellipse), size = 0.6, alpha = 0.7) +
  scale_color_manual(values = c("FALSE" = "grey80", "TRUE" = "red")) +
  geom_path(data = ellipse_outline, aes(x = x, y = y), inherit.aes = FALSE,
            color = "blue", linewidth = 0.8) +
  labs(title = paste0("Rotated oval gate (", angle_deg, "\u00b0 CCW) -- adjust as needed"),
       x = "nclTEC module score", y = "PolyKRT module score", color = "In gate") +
  theme_minimal(base_size = 11)
print(p_ellipse_check)
ggsave(paste0(DATE_PREFIX, "-central-cloud-ellipse-gate-rotated-check.pdf"),
       p_ellipse_check, width = 8, height = 7)

samples_central_cloud <- subset(samples, cells = scatter_df_ellipse$cell[scatter_df_ellipse$in_ellipse])
cat("\nSubsetted Seurat object cell count:", ncol(samples_central_cloud), "\n")



# -----------------------------------------------------------------------------
# Subset cells in BOTH the rotated ellipse gate AND annotated as BMC TEC
# (self-contained -- re-derives the ellipse gate and BMC TEC filter from
# scratch, doesn't assume prior objects are still in memory)
# -----------------------------------------------------------------------------
center_x <- 1.5
center_y <- 0.4
radius_x <- 0.8
radius_y <- 0.9
angle_deg <- 25
angle_rad <- angle_deg * pi / 180

scatter_df_ellipse <- scatter_df %>%
  tibble::rownames_to_column("cell") %>%
  mutate(
    x_centered = nclTEC_score_2d - center_x,
    y_centered = PolyKRT_score_2d - center_y,
    x_rot = x_centered * cos(-angle_rad) - y_centered * sin(-angle_rad),
    y_rot = x_centered * sin(-angle_rad) + y_centered * cos(-angle_rad),
    ellipse_dist = (x_rot / radius_x)^2 + (y_rot / radius_y)^2,
    in_ellipse = ellipse_dist <= 1
  )

# sanity check: confirm ellipse gate alone still captures a non-trivial subset
n_in_ellipse <- sum(scatter_df_ellipse$in_ellipse)
cat("Cells within ellipse gate (before BMC TEC filter):", n_in_ellipse, "\n")
stopifnot(n_in_ellipse > 0)

ellipse_gated_cells <- scatter_df_ellipse$cell[scatter_df_ellipse$in_ellipse]

# -----------------------------------------------------------------------------
# Intersect with BMC TEC annotation
# -----------------------------------------------------------------------------
bmc_tec_cells <- colnames(samples)[samples$cell_type_unintegrated == "BMC TEC"]

# sanity check: confirm BMC TEC cells exist in the object
stopifnot(length(bmc_tec_cells) > 0)
cat("Total BMC TEC cells in dataset:", length(bmc_tec_cells), "\n")

ellipse_and_bmc_cells <- intersect(ellipse_gated_cells, bmc_tec_cells)

# sanity check: confirm the intersection is non-empty and smaller than either input set
cat("Cells in BOTH ellipse gate AND BMC TEC:", length(ellipse_and_bmc_cells), "\n")
stopifnot(length(ellipse_and_bmc_cells) > 0)
stopifnot(length(ellipse_and_bmc_cells) <= min(length(ellipse_gated_cells), length(bmc_tec_cells)))

# -----------------------------------------------------------------------------
# Visual check: highlight this intersection specifically
# -----------------------------------------------------------------------------
theta <- seq(0, 2 * pi, length.out = 200)
ellipse_outline <- data.frame(
  x_unrot = radius_x * cos(theta), y_unrot = radius_y * sin(theta)
) %>%
  mutate(x = center_x + x_unrot * cos(angle_rad) - y_unrot * sin(angle_rad),
         y = center_y + x_unrot * sin(angle_rad) + y_unrot * cos(angle_rad))

scatter_df_check <- scatter_df %>%
  tibble::rownames_to_column("cell") %>%
  mutate(in_both = cell %in% ellipse_and_bmc_cells)

p_both_check <- ggplot(scatter_df_check, aes(x = nclTEC_score_2d, y = PolyKRT_score_2d)) +
  geom_point(aes(color = in_both), size = 0.6, alpha = 0.7) +
  scale_color_manual(values = c("FALSE" = "grey80", "TRUE" = "red")) +
  geom_path(data = ellipse_outline, aes(x = x, y = y), inherit.aes = FALSE,
            color = "blue", linewidth = 0.8) +
  labs(title = "Cells in ellipse gate AND annotated as BMC TEC",
       x = "nclTEC module score", y = "PolyKRT module score", color = "In both") +
  theme_minimal(base_size = 11)
print(p_both_check)
ggsave(paste0(DATE_PREFIX, "-ellipse-AND-BMCTEC-check.pdf"), p_both_check, width = 8, height = 7)

# -----------------------------------------------------------------------------
# Final Seurat subset
# -----------------------------------------------------------------------------
samples_ellipse_bmc <- subset(samples, cells = ellipse_and_bmc_cells)
cat("\nFinal subsetted Seurat object cell count:", ncol(samples_ellipse_bmc), "\n")
stopifnot(ncol(samples_ellipse_bmc) == length(ellipse_and_bmc_cells))



# -----------------------------------------------------------------------------
# Annotate ellipse-gated BMC TEC cells directly in samples@meta.data
# -----------------------------------------------------------------------------
center_x <- 1.5; center_y <- 0.4; radius_x <- 0.8; radius_y <- 0.9
angle_rad <- 25 * pi / 180

meta_df <- FetchData(samples, vars = c("nclTEC_score_2d", "PolyKRT_score_2d", "cell_type_unintegrated")) %>%
  tibble::rownames_to_column("cell") %>%
  mutate(
    x_rot = (nclTEC_score_2d - center_x) * cos(-angle_rad) - (PolyKRT_score_2d - center_y) * sin(-angle_rad),
    y_rot = (nclTEC_score_2d - center_x) * sin(-angle_rad) + (PolyKRT_score_2d - center_y) * cos(-angle_rad),
    in_ellipse = (x_rot / radius_x)^2 + (y_rot / radius_y)^2 <= 1,
    in_group = in_ellipse & cell_type_unintegrated == "BMC TEC"
  )

samples$ellipse_bmc_group <- ifelse(meta_df$in_group, "Ellipse_BMC_TEC", "Other")

# sanity check: counts sum to total, matches intersection size
stopifnot(sum(table(samples$ellipse_bmc_group)) == ncol(samples))
print(table(samples$ellipse_bmc_group))

# sanity check: plot the annonated "Ellipse BMC TEC" on the module score scatter
# plot to confirm its the correct cells being annoteated as ellipse bmc TEC
p_verify <- ggplot(meta_df, aes(x = nclTEC_score_2d, y = PolyKRT_score_2d)) +
  geom_point(aes(color = in_group), size = 0.6, alpha = 0.7) +
  scale_color_manual(values = c("FALSE" = "grey80", "TRUE" = "red")) +
  labs(title = "Sanity check: Ellipse_BMC_TEC cells (red) on module score scatter",
       x = "nclTEC module score", y = "PolyKRT module score", color = "In group")
print(p_verify)
ggsave(paste0(DATE_PREFIX, "-ellipse-bmc-sanity-check.pdf"), p_verify, width = 8, height = 7)


# -----------------------------------------------------------------------------
# dimplot highlighting ellipse bmc tec
# -----------------------------------------------------------------------------
p_dimplot_check <- DimPlot(
  samples, reduction = "umap_unintegrated",
  cells.highlight = colnames(samples)[samples$ellipse_bmc_group == "Ellipse_BMC_TEC"],
  sizes.highlight = 0.6, pt.size = 0.3
) +
  scale_color_manual(values = c("grey85", "red"), labels = c("Other", "Ellipse_BMC_TEC")) +
  labs(title = "Ellipse-gated BMC TEC cells on whole-dataset UMAP")
print(p_dimplot_check)
ggsave(paste0(DATE_PREFIX, "-ellipse-bmc-dimplot-check.pdf"), p_dimplot_check, width = 7, height = 6)


# -----------------------------------------------------------------------------
# dimplot highlighting ellipse bmc tec
# -----------------------------------------------------------------------------
cross_tab_sample <- table(samples$ellipse_bmc_group, samples$sample.id)
print(cross_tab_sample)
print(round(100 * prop.table(cross_tab_sample, margin = 2), 1))

