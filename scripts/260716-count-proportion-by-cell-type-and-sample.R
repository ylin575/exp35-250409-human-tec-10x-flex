# =============================================================================
# Cell type composition of the sorted CD205+ TEC sample
# Answers Reviewer 2's comment about mTEC I appearing at the bottom of the
# cTEC cluster in the CD205+ sample.
#
# Sample structure (verified from meta.data):
#   ht67-cd205pos  <- sorted CD205+ (target of the reviewer question)
#   ht67-cd205neg  <- sorted CD205- (same donor, opposite gate -- best control)
#   ht70, ht71     <- unsorted "All TECs" (different donors)
# =============================================================================

library(Seurat)
library(dplyr)
library(tidyr)

tec <- readRDS("data/rds-objects/260501-after-annotation-1.rds")

# -----------------------------------------------------------------------------
# 1. Composition of CD205+ sample specifically -- the direct answer to Reviewer 2
# -----------------------------------------------------------------------------
cd205pos_composition <- tec@meta.data %>%
  filter(sample.id == "ht67-cd205pos") %>%
  count(cell_type, name = "n_cells") %>%
  mutate(percent = round(100 * n_cells / sum(n_cells), 2)) %>%
  arrange(desc(n_cells))

cd205pos_composition <- bind_rows(
  cd205pos_composition,
  tibble(cell_type = "TOTAL",
         n_cells   = sum(cd205pos_composition$n_cells),
         percent   = 100)
)

cat("=== CD205+ (ht67-cd205pos) composition ===\n")
print(cd205pos_composition)

write.csv(cd205pos_composition,
          "CD205pos_composition.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 2. All four samples side-by-side -- context for the response letter
# -----------------------------------------------------------------------------
composition_all <- tec@meta.data %>%
  count(sample.id, cell_type, name = "n_cells") %>%
  group_by(sample.id) %>%
  mutate(percent = round(100 * n_cells / sum(n_cells), 2)) %>%
  ungroup()

# Long format for the CSV
write.csv(composition_all,
          "composition_all_samples_long.csv", row.names = FALSE)

# Wide format for at-a-glance comparison (percent per sample per cell type)
composition_wide <- composition_all %>%
  select(sample.id, cell_type, percent) %>%
  pivot_wider(names_from = sample.id, values_from = percent, values_fill = 0) %>%
  select(cell_type,
         `ht67-cd205pos`, `ht67-cd205neg`, ht70, ht71)   # column order

cat("\n=== Percent composition, all four samples (columns) ===\n")
print(composition_wide)

write.csv(composition_wide,
          "composition_all_samples_wide.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 3. The paired comparison Reviewer 2's argument turns on:
#    mTEC I proportion in CD205+ vs CD205- from the SAME donor
# -----------------------------------------------------------------------------
paired_ht67 <- composition_wide %>%
  select(cell_type, `ht67-cd205pos`, `ht67-cd205neg`) %>%
  mutate(
    diff_pos_minus_neg = `ht67-cd205pos` - `ht67-cd205neg`
  )

cat("\n=== HT67 paired (CD205+ vs CD205-) composition ===\n")
print(paired_ht67)

write.csv(paired_ht67, "ht67_paired_composition.csv", row.names = FALSE)