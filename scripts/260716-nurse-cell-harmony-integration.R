# =============================================================================
# Nurse cell subclustering WITH HARMONY INTEGRATION
#
# Previous attempt without integration showed near-total correspondence between
# subcluster identity and donor of origin (ht67-cd205pos, ht70, ht71 each
# dominating a different cluster). This script integrates over sample.id via
# Harmony to remove that donor structure, then re-clusters on the corrected
# embedding to test whether biology-driven subpopulations exist underneath.
# =============================================================================

# install.packages("harmony")   # if not already
library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(patchwork)
library(scCustomize)

# Date prefix for all output filenames (yymmdd-)
DATE_PREFIX <- format(Sys.time(), "%y%m%d-%H%M")

tec <- readRDS("data/rds-objects/260501-after-annotation-1.rds")

sample_metadata <- data.frame(
  sample.id   = c("ht67-cd205pos", "ht67-cd205neg", "ht70", "ht71"),
  age_months  = c(3, 3, 12, 48),
  sex         = c("F", "F", "M", "F"),
  institution = c("Mount Sinai", "Mount Sinai", "Mount Sinai", "Northwell"),
  stringsAsFactors = FALSE
)

tec@meta.data <- tec@meta.data %>%
  left_join(sample_metadata, by = "sample.id")
rownames(tec@meta.data) <- colnames(tec)   # restore rownames -- dplyr strips them

# Verify
table(tec$sample.id, tec$age_months, tec$sex, tec$institution)

# -----------------------------------------------------------------------------
# 1. Subset to nurse cells, wipe stale parent layers/reductions
# -----------------------------------------------------------------------------
nurse <- subset(tec, subset = cell_type == "Nurse")

cat("Total nurse cells:", ncol(nurse), "\n")
cat("Per-sample distribution:\n")
print(table(nurse$sample.id))

DefaultAssay(nurse) <- "RNA"
nurse <- DietSeurat(
  nurse,
  layers    = c("counts", "data"),
  dimreducs = NULL,
  graphs    = NULL,
  misc      = FALSE
)

# -----------------------------------------------------------------------------
# 2. Standard preprocessing through PCA
# -----------------------------------------------------------------------------
nurse <- nurse %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(nfeatures = 2000, verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE)

# -----------------------------------------------------------------------------
# 3. HARMONY INTEGRATION -- correct on sample.id
# -----------------------------------------------------------------------------
# Harmony iteratively adjusts PCA embeddings so cells cluster by biology
# rather than by batch/donor. Output is a new reduction called "harmony"
# with the same dimensionality as the input PCA; use it in place of "pca"
# for all downstream steps.
#
# group.by.vars = "sample.id" tells Harmony which variable(s) to correct.
# theta controls correction strength (default 2). Higher = more aggressive
# correction, but risks removing real biology if it correlates with batch.

nurse <- RunHarmony(
  nurse,
  group.by.vars = "sample.id",
  reduction     = "pca",
  dims.use      = 1:30,
  reduction.save = "harmony",
  verbose       = FALSE
)

# -----------------------------------------------------------------------------
# 4. Downstream steps ON THE HARMONY EMBEDDING (not PCA)
# -----------------------------------------------------------------------------
N_DIMS <- 20   # same rationale as before; harmony preserves original dim count

nurse <- nurse %>%
  FindNeighbors(reduction = "harmony", dims = 1:N_DIMS, verbose = FALSE) %>%
  RunUMAP(reduction = "harmony", dims = 1:N_DIMS, verbose = FALSE)

# -----------------------------------------------------------------------------
# 5. Sanity check FIRST -- did Harmony actually mix the samples?
# -----------------------------------------------------------------------------
# If sample.id is now well-distributed across the UMAP, integration worked.
# If samples still segregate visibly, either theta needs to increase or the
# donor effect is so strong it can't be corrected without losing biology.

p_mix_check <- DimPlot_scCustom(nurse, group.by = "sample.id") +
  labs(title = "After Harmony -- samples should now be mixed") +
  theme(legend.position = "right")

print(p_mix_check)
ggsave(paste0(DATE_PREFIX, "-nurse-harmony-sample-mixing-check.pdf"), p_mix_check,
       width = 6, height = 5)

# -----------------------------------------------------------------------------
# 6. Cluster at multiple resolutions on the integrated embedding
# -----------------------------------------------------------------------------
resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0)
for (r in resolutions) {
  nurse <- FindClusters(nurse, resolution = r, verbose = FALSE)
}

cat("\n=== Cluster counts at each resolution (post-Harmony) ===\n")
for (r in resolutions) {
  col <- paste0("RNA_snn_res.", r)
  cat(sprintf("res %.1f  ->  %d clusters\n", r, length(unique(nurse[[col]][[1]]))))
}

p_res <- wrap_plots(
  lapply(resolutions, function(r) {
    col <- paste0("RNA_snn_res.", r)
    DimPlot_scCustom(nurse, group.by = col, label = TRUE, repel = TRUE) +
      labs(title = paste0("res = ", r)) +
      NoLegend()
  }),
  ncol = 3
)
print(p_res)
ggsave(paste0(DATE_PREFIX, "-nurse-harmony-subclusters-resolutions.pdf"), p_res,
       width = 12, height = 8)

# Choose resolution -- pick the lowest that shows visually distinct groups
CHOSEN_RES <- 0.4                                 # <-- adjust after inspecting
Idents(nurse) <- paste0("RNA_snn_res.", CHOSEN_RES)

# -----------------------------------------------------------------------------
# 7. Re-run the donor cross-tab -- confirm integration worked at cluster level
# -----------------------------------------------------------------------------
cat("\n=== Post-Harmony subcluster vs. sample.id ===\n")
print(table(Idents(nurse), nurse$sample.id))

# INTERPRETATION:
#   Each cluster should now contain cells from MULTIPLE samples in roughly
#   representative proportions. If clusters are still >80% single-sample,
#   integration didn't fully resolve the donor effect for those clusters --
#   they may represent genuine donor-specific biology or residual batch.

p_check <- (
  DimPlot_scCustom(nurse, group.by = paste0("RNA_snn_res.", CHOSEN_RES),
                   label = TRUE) + labs(title = "Post-Harmony subclusters")
) | (
  DimPlot_scCustom(nurse, group.by = "sample.id") + labs(title = "By donor/sort")
)
print(p_check)
ggsave(paste0(DATE_PREFIX, "-nurse-harmony-subclusters-vs-donor.pdf"), p_check,
       width = 12, height = 5)








# =============================================================================
# SECTION 2: Age-coloring diagnostic on the nurse UMAP
# =============================================================================

# Donor age lookup -- numeric (months) for a continuous, ordered color scale
age_lookup <- c(
  "ht67-cd205pos" = 3,    # 3 months
  "ht67-cd205neg" = 3,    # 3 months (same donor)
  "ht70"          = 12,   # 1 year
  "ht71"          = 48    # 4 years
)

nurse$age_months <- unname(age_lookup[nurse$sample.id])

table(nurse$sample.id, nurse$age_months, useNA = "ifany")

# Sanity check: confirm every cell got an age assigned (no NAs)
cat("\n=== Age assignment check ===\n")
print(table(nurse$sample.id, nurse$age_months, useNA = "ifany"))

# Plot: nurse UMAP colored by age (continuous scale) -- uses whatever UMAP
# reduction currently exists on `nurse` (pre-Harmony PCA-based UMAP if this
# runs before Section 3, or update reduction= if you've already run Harmony)
p_age <- FeaturePlot_scCustom(
  seurat_object = nurse,
  reduction     = "umap",
  features      = "age_months",
  colors_use    = viridis_light_high,
  order         = TRUE,
  pt.size       = 0.6
) +
  labs(title = "Nurse UMAP colored by donor age (months)")

print(p_age)
ggsave(paste0(DATE_PREFIX, "-nurse-umap-by-age.pdf"), p_age,
       width = 6, height = 5)

# Side-by-side: age vs. sample.id, for direct visual comparison
p_age_vs_sample <- p_age | (
  DimPlot_scCustom(nurse, group.by = "sample.id") +
    labs(title = "By sample.id")
)
print(p_age_vs_sample)
ggsave(paste0(DATE_PREFIX, "-nurse-umap-age-vs-sample.pdf"), p_age_vs_sample,
       width = 12, height = 5)

# INTERPRETATION:
#   A clean, monotonic age gradient across the UMAP (e.g. one pole = 3mo,
#   opposite pole = 48mo, with 12mo in between) supports age as a real
#   biological driver of structure -- argues AGAINST aggressively increasing
#   Harmony theta to force further mixing, since you'd be erasing age biology.
#   If age coloring looks just as scattered/non-gradient as sample.id, age
#   alone doesn't explain the structure -- individual donor idiosyncrasy or
#   residual technical batch remain live explanations.


# =============================================================================
# SECTION 3: Harmony theta sweep with quantitative mixing metric
# =============================================================================
# Starts from `nurse` AFTER RunPCA but BEFORE any RunHarmony call.
# Loops over theta values, reruns Harmony + neighbors + UMAP + clustering at
# a FIXED resolution for each, and computes a chi-square-based skew metric
# per theta (no extra packages required).
#
# OPTIONAL UPGRADE: install.packages("lisi") [github.com/immunogenomics/LISI]
# for the field-standard iLISI mixing score instead of/alongside this metric.

thetas <- c(0, 1, 2, 4, 6, 8)
N_DIMS <- 20                      # match what you used before
FIXED_RES <- 0.2                  # match your chosen nurse resolution

baseline_props <- prop.table(table(nurse$sample.id))  # for chi-square expected

theta_results <- list()
theta_plots <- list()

for (th in thetas) {
  cat("\n--- theta =", th, "---\n")
  
  red_name <- paste0("harmony_theta", th)
  nurse <- RunHarmony(
    nurse,
    group.by.vars  = "sample.id",
    reduction      = "pca",
    dims.use       = 1:30,
    theta          = th,
    reduction.save = red_name,
    verbose        = FALSE
  )
  
  nurse <- FindNeighbors(nurse, reduction = red_name, dims = 1:N_DIMS,
                         graph.name = paste0("nn_theta", th), verbose = FALSE)
  nurse <- FindClusters(nurse, resolution = FIXED_RES,
                        graph.name = paste0("nn_theta", th),
                        cluster.name = paste0("clusters_theta", th),
                        verbose = FALSE)
  nurse <- RunUMAP(nurse, reduction = red_name, dims = 1:N_DIMS,
                   reduction.name = paste0("umap_theta", th), verbose = FALSE)
  
  # --- Quantitative mixing metric: per-cluster chi-square vs baseline ---
  clust_col <- paste0("clusters_theta", th)
  tab <- table(nurse[[clust_col]][[1]], nurse$sample.id)
  
  chisq_per_cluster <- apply(tab, 1, function(row) {
    expected <- baseline_props * sum(row)
    suppressWarnings(chisq.test(x = row, p = baseline_props)$p.value)
  })
  
  pct_clusters_skewed <- round(100 * mean(chisq_per_cluster < 0.05), 1)
  n_clusters <- nrow(tab)
  
  theta_results[[as.character(th)]] <- tibble(
    theta = th,
    n_clusters = n_clusters,
    pct_clusters_significantly_skewed = pct_clusters_skewed
  )
  
  # Visual: UMAP colored by sample.id at this theta
  theta_plots[[as.character(th)]] <- DimPlot_scCustom(
    nurse, reduction = paste0("umap_theta", th), group.by = "sample.id"
  ) +
    labs(title = paste0("theta = ", th)) +
    NoLegend() +
    theme(plot.title = element_text(size = 10))
}

# -----------------------------------------------------------------------------
# Summary table across all thetas
# -----------------------------------------------------------------------------
theta_summary <- bind_rows(theta_results)

cat("\n=== Theta sweep summary ===\n")
cat("pct_clusters_significantly_skewed = % of clusters whose sample\n")
cat("composition significantly differs from the dataset-wide baseline\n")
cat("(chi-square p < 0.05). LOWER = better mixing.\n\n")
print(theta_summary)

write.csv(theta_summary, paste0(DATE_PREFIX, "-harmony-theta-sweep-summary.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# Visual grid across all thetas
# -----------------------------------------------------------------------------
p_theta_grid <- wrap_plots(theta_plots, ncol = 3) +
  plot_layout(guides = "collect")

print(p_theta_grid)
ggsave(paste0(DATE_PREFIX, "-harmony-theta-sweep-umaps.pdf"), p_theta_grid,
       width = 12, height = 8)

# -----------------------------------------------------------------------------
# INTERPRETATION
# -----------------------------------------------------------------------------
cat("\n", strrep("=", 70), "\n")
cat("HOW TO READ THE SWEEP:\n\n")
cat("- If pct_clusters_skewed drops sharply and plateaus at some theta,\n")
cat("  that plateau is a reasonable choice -- pushing theta higher past\n")
cat("  the plateau adds correction strength without improving mixing,\n")
cat("  which risks over-correction (erasing real biology, e.g. age).\n\n")
cat("- If pct_clusters_skewed never drops much regardless of theta, the\n")
cat("  residual structure is likely NOT a correctable technical batch\n")
cat("  effect -- consistent with the age-confound explanation from\n")
cat("  Section 2, and argues against chasing higher theta further.\n\n")
cat("- Cross-reference against Section 2's age plot: does the theta value\n")
cat("  that best 'mixes' sample.id also blur the age gradient? If so,\n")
cat("  that theta is removing biology you may want to keep.\n")
cat(strrep("=", 70), "\n")

# Save the object with all theta variants retained for further inspection
saveRDS(nurse, paste0(DATE_PREFIX, "-nurse-theta-sweep.rds"))





















# -----------------------------------------------------------------------------
# 8. Find markers on integrated subclusters
# -----------------------------------------------------------------------------
# NOTE: FindAllMarkers uses the RNA data layer (uncorrected), NOT the Harmony
# embedding -- this is standard practice. Harmony only corrects the reduced-
# dimension embedding used for neighborhood detection; markers should still
# come from the actual expression data.

nurse_markers <- FindAllMarkers(
  nurse,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25,
  verbose         = FALSE
)

top_markers <- nurse_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 15) %>%
  ungroup()

write.csv(nurse_markers, paste0(DATE_PREFIX, "-nurse-harmony-markers-ALL.csv"), row.names = FALSE)
write.csv(top_markers,   paste0(DATE_PREFIX, "-nurse-harmony-markers-TOP15.csv"), row.names = FALSE)

cat("\n=== Top markers per subcluster (post-Harmony) ===\n")
print(top_markers %>% select(cluster, gene, avg_log2FC, pct.1, pct.2, p_val_adj))

# -----------------------------------------------------------------------------
# 9. Noise diagnostic
# -----------------------------------------------------------------------------
noise_patterns <- c("^RPS", "^RPL", "^MT-", "^HSP", "^HIST",
                    "^MKI67$", "^TOP2A$", "^CDK1$", "^MCM", "^CCN")
noise_regex <- paste(noise_patterns, collapse = "|")

diagnostic <- top_markers %>%
  mutate(is_noise = grepl(noise_regex, gene)) %>%
  group_by(cluster) %>%
  summarise(
    n_top       = n(),
    n_noise     = sum(is_noise),
    pct_noise   = round(100 * n_noise / n_top, 1),
    noise_genes = paste(gene[is_noise], collapse = ", "),
    .groups = "drop"
  )

cat("\n=== Noise diagnostic (post-Harmony) ===\n")
print(diagnostic)

# -----------------------------------------------------------------------------
# 10. Marker dotplot
# -----------------------------------------------------------------------------
top5 <- top_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5) %>%
  pull(gene) %>%
  unique()

p_dot <- DotPlot_scCustom(nurse, features = top5,
                          group.by = paste0("RNA_snn_res.", CHOSEN_RES)) +
  RotatedAxis() +
  labs(title = "Top-5 markers per nurse subcluster (post-Harmony)")
print(p_dot)
ggsave(paste0(DATE_PREFIX, "-nurse-harmony-marker-dotplot.pdf"), p_dot,
       width = max(8, length(top5) * 0.35), height = 5)

# -----------------------------------------------------------------------------
# 11. Save
# -----------------------------------------------------------------------------
saveRDS(nurse, paste0(DATE_PREFIX, "-nurse-harmony-subclustered.rds"))

cat("\n=== Outputs (prefixed with", DATE_PREFIX, ") ===\n")
cat("  -nurse-harmony-sample-mixing-check.pdf     -- did Harmony work?\n")
cat("  -nurse-harmony-subclusters-resolutions.pdf -- pick a resolution\n")
cat("  -nurse-harmony-subclusters-vs-donor.pdf    -- integration at cluster level\n")
cat("  -nurse-harmony-markers-ALL.csv             -- full marker table\n")
cat("  -nurse-harmony-markers-TOP15.csv           -- top 15 per subcluster\n")
cat("  -nurse-harmony-marker-dotplot.pdf          -- visual marker summary\n")
cat("  -nurse-harmony-subclustered.rds            -- integrated object\n")