# load libraries

library(MuSiC) # for deconvolution analysis
library(readr) # for importing tsv file into data frame
library(Biobase) # contains the ExpressionSet() function for bulk rnaseq data
library(SingleCellExperiment) # for loading scrnaseq dataset
library(Seurat) # 
library(tidyverse) #
library(biomaRt) # needed to convert gene id to gene name for the bulk dataset
library(cowplot) # to plot the jitter plots

# load bulk rnaseq raw count into ExpressionSet object

# read in tsv count matrix into data frame
bulk_counts <- read_tsv("cd49f_ctec_analysis/counts/all.tsv")
# add column of sum ctec and mtec gene counts
bulk_counts <- bulk_counts %>%
  # mean across replicates
  mutate(ctec_mean = (ctec1 + ctec2 + ctec3)/3,
         mtec_mean = (mtec1 + mtec2 + mtec3)/3,
         # remove version number (.xxx) of gene id
         gene_id=gsub("\\..*", "", gene)) %>%
  # arrange gene column
  arrange(gene)


#bulk_counts <- bulk_counts %>%
#  mutate(mtec_sum_across = across(starts_with("mtec"), ~ sum(., na.rm=TRUE)),
#        mtec_sum = mtec1 + mtec2 + mtec3)
# remove the ".xxx" from the gene names
#bulk_counts$gene <- gsub("\\..*", "", bulk_counts$gene)
# arrange gene column
#bulk_counts <- arrange(bulk_counts, gene)


# need to convert gene id to gene name to be compatible with the scrnaseq reference dataset
# download human gene dataset from ensembl using biomart

ensembl <- useMart("ensembl")
dataset <- useDataset("hsapiens_gene_ensembl", mart = ensembl)

# Query for gene names
gene_info <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = bulk_counts$gene_id,
  mart = dataset
)
# save gene info
saveRDS(gene_info, "gene_info.rds")
gene_info <- readRDS("gene_info.rds")

# Print result
print(gene_info)
# sort gene info alphabetically
gene_info_sorted <- gene_info |> arrange(ensembl_gene_id)

# filter bulk_counts to those only in gene_info
bulk_counts <- bulk_counts[bulk_counts$gene_id %in% gene_info[,1], ]

## format duplicate gene ids
# do a sanity check to find out what duplicate ids are
unique_duplicates <- unique(bulk_counts$gene_id[duplicated(bulk_counts$gene_id)])
# seems like duplicates are on pseduo-autosomal region; let's get the mean
bulk_counts |> filter(gene_id %in% unique_duplicates)

bulk_counts <- bulk_counts |>
  # group by unique gene_id
  group_by(gene_id) |>
  # how do we want the rows that belong to unique gene_ids to be summarised:
  # find all columns that have doubles (ie all the ones with counts)
  # for gene_ids with multiple entries, get the mean for each column
  summarise(across(where(is.double), ~ mean(.x, na.rm = TRUE)))

# extract row entries of the duplicate genes from dataset
#bulk_counts_duplicate_genes <- bulk_counts[bulk_counts$gene %in% unique_duplicates,]

# remove row entries that have duplicate entries
#bulk_counts <- bulk_counts[!bulk_counts$gene %in% unique_duplicates, ]

# remove same genes from gene_info_sorted so the row numbers match for cbind
# gene_info_sorted <- gene_info_sorted[!gene_info_sorted$ensembl_gene_id %in% unique_duplicates, ]

# combine gene_info_sorted and bulk_counts
a <- cbind(gene_info_sorted, bulk_counts)
# check identical column 1 and 3 of a
identical(a[,1], a[,3])

# select only the columns needed for MuSiC
bulk_avgs <- a[,c(2,10,11)]
bulk_indiv <- a[,c(2,4,5,6,7,8,9)]

# find duplicate gene names to remove before converting gene name to rowname
duplicates_gene_info <- unique(gene_info$external_gene_name[duplicated(gene_info$external_gene_name)])
  # looks like there are 19256 rows that have duplicate gene names
  # unique = 79

# remove entries with duplicate gene names
bulk_avgs <- bulk_avgs[!bulk_avgs$external_gene_name %in% duplicates_gene_info,]
bulk_indiv <- bulk_indiv[!bulk_indiv$external_gene_name %in% duplicates_gene_info,]

# add genes as rownames
rownames(bulk_avgs) <- bulk_avgs$external_gene_name
rownames(bulk_indiv) <- bulk_indiv$external_gene_name

# remove external_gene_name column before converting to matrix
bulk_avgs <- bulk_avgs[,c(2,3)]
bulk_indiv <- bulk_indiv[,c(2,3,4,5,6,7)]

# need to convert to matrix before using ExpressionSet()
bulk_avgs <- as.matrix(bulk_avgs)
bulk_avgs <- ExpressionSet(assayData = bulk_avgs)
bulk_indiv <- as.matrix(bulk_indiv)
bulk_indiv <- ExpressionSet(assayData = bulk_indiv)


# load single cell data: will use parse data with labels transferred from multiome
hthy2wi <- readRDS("240205_all_7samples_integrated-pre_label_transfer.rds")

# perform label transfer from multiome to parse dataset
tec.ref <- subset(hthy2wi, orig.ident == "SeuratProject") # the multiome dataset
tec.query <- subset(hthy2wi, orig.ident == "hthy") # the parse dataset

tec.anchors <- FindTransferAnchors(reference = tec.ref, 
                                   query = tec.query,
                                   dims = 1:30, 
                                   reference.reduction = "pca")
predictions <- TransferData(anchorset = tec.anchors, 
                            refdata = tec.ref$celltype, 
                            dims = 1:30)
tec.query <- AddMetaData(tec.query, metadata = predictions)

# Unimodal UMAP Projection and MapQuery
tec.ref <- RunUMAP(tec.ref, dims = 1:30, 
                   reduction = "integrated.cca", 
                   return.model = TRUE)
tec.query <- MapQuery(anchorset = tec.anchors, 
                      reference = tec.ref, 
                      query = tec.query,
                      refdata = list(celltype = "celltype"), 
                      reference.reduction = "pca", reduction.model = "umap")

# subset parse to get only ht12 and ht14
t <- subset(tec.query, sample == c("Ht12", "Ht14"))
# convert parse dataset to SingleCellExperiment object
tec.query.sce <- as.SingleCellExperiment(tec.query) # parse dataset
tec.ref.sce <- as.SingleCellExperiment(tec.ref) # multiome dataset
t.sce <- as.SingleCellExperiment(t) # parse dataset: ht12 and ht14 only


# extract cell types from single cell datasets
celltype.parse <- unique(tec.query@meta.data$predicted.celltype)
celltype.multiome <- unique(tec.ref@meta.data$celltype)


# use MuSiC

# extract expression data from bulk ExpressionSet object
bulk_avgs.mtx = exprs(bulk_avgs)
bulk_indiv.mtx = exprs(bulk_indiv)

# run MuSiC using parse dataset (tec.query.sce) as reference
est.prop.avgs = music_prop(bulk.mtx = bulk_avgs.mtx, sc.sce = tec.query.sce,
                          clusters = 'predicted.celltype',
                      samples = 'samples', select.ct = celltype.parse)
est.prop.indiv = music_prop(bulk.mtx = bulk_indiv.mtx, sc.sce = tec.query.sce,
                            clusters = 'predicted.celltype',
                          samples = 'samples', select.ct = celltype.parse)
head(est.prop.avgs)
head(est.prop.indiv)

# run MuSiC using parse dataset: ht12 and ht14 only (t.sce) as reference
est.prop.avgs = music_prop(bulk.mtx = bulk_avgs.mtx, sc.sce = t.sce,
                          clusters = 'predicted.celltype',
                          samples = 'samples', select.ct = celltype.parse)
est.prop.indiv = music_prop(bulk.mtx = bulk_indiv.mtx, sc.sce = t.sce,
                            clusters = 'predicted.celltype',
                            samples = 'samples', select.ct = celltype.parse)
head(est.prop.avgs)
head(est.prop.indiv)

# run MuSiC using multiome dataset (tec.ref.sce) as reference
est.prop.avgs = music_prop(bulk.mtx = bulk_avgs.mtx, sc.sce = tec.ref.sce, clusters = 'celltype',
                          samples = 'samples', select.ct = celltype.multiome)
est.prop.indiv = music_prop(bulk.mtx = bulk_indiv.mtx, sc.sce = tec.ref.sce, clusters = 'celltype',
                             samples = 'samples', select.ct = celltype.multiome)
head(est.prop.avgs)
head(est.prop.indiv)


# plotting jitter plots

jitter.fig.avgs = Jitter_Est(data.matrix(est.prop.sum$Est.prop.weighted),
                        method.name = 'MuSiC', title = 'Jitter plot of Est Proportions')

jitter.fig.indiv = Jitter_Est(data.matrix(est.prop.indiv$Est.prop.weighted),
                            method.name = 'MuSiC', title = 'Jitter plot of Est Proportions')

plot_grid(jitter.fig.sum)
plot_grid(jitter.fig.indiv)



# different plotting, using tec.ref.sce as reference

# output data
df <- data.frame(
  row.names = c("ctec1", "ctec2", "ctec3", "mtec1", "mtec2", "mtec3"),
  mTEC_I = c(0.0000000, 0.0000000, 0.0000000, 0.2534117, 0.3080576, 0.3743840),
  neuro = c(0.0000000, 0.0000000, 0.0000000, 0.0142836, 0.0000000, 0.0000000),
  neuroendocrine_I = c(0.00000000, 0.01698763, 0.00000000, 0.08591833, 0.05350527, 0.10251167),
  mcTEC = c(0.5941152, 0.5472839, 0.4688312, 0.0000000, 0.0000000, 0.0000000),
  sTEC_II = c(0.05686663, 0.05922092, 0.07237718, 0.06410879, 0.06894198, 0.07796896),
  muscle = c(0.0334431140, 0.0001394828, 0.0041080313, 0.0040013988, 0.0295649559, 0.0377315275),
  neuroendocrine_II = c(0.0001368444, 0.0005560554, 0.0060538083, 0.0537450774, 0.0000000000, 0.0000000000),
  ionocyte_II = c(0.03679028, 0.07310675, 0.06819035, 0.03063144, 0.08203186, 0.08934798),
  ionocyte_I = c(0, 0, 0, 0, 0, 0),
  mTEC_II = c(0, 0, 0, 0, 0, 0),
  tuft = c(0, 0, 0, 0, 0, 0),
  sTEC_I = c(0.0345251381, 0.0100835103, 0.0612530555, 0.0000000000, 0.0007995067, 0.0030092502),
  transit_amplifying = c(0.1328722, 0.1657999, 0.1730291, 0.3978898, 0.2272830, 0.1504242),
  cTEC_late = c(0.07582003, 0.04097879, 0.01625137, 0.00000000, 0.00000000, 0.00000000),
  cTEC_early = c(0.00000000, 0.05351126, 0.05432203, 0.00000000, 0.02158990, 0.01231656),
  heteroTEC = c(0.03543052, 0.03233180, 0.07558385, 0.09600995, 0.20822590, 0.15230584)
)


# Convert to long format
df_long <- df %>%
  tibble::rownames_to_column(var = "Sample") %>%
  pivot_longer(cols = -Sample, names_to = "Type", values_to = "Value")


# Calculate the common y-axis limits
y_limits <- range(df_long$Value, na.rm = TRUE)


create_plot <- function(data, column_name) {
  ggplot(data, aes(x = Sample, y = Value, label = Sample)) +
    geom_point() +
    geom_text(vjust = -0.5) +
    ylim(y_limits) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = column_name, x = "Sample", y = "Value")
}



# Generate a list of plots
plot_list <- lapply(names(df)[-1], function(col) {
  df_subset <- df_long %>% filter(Type == col)
  create_plot(df_subset, col)
})

# Using cowplot
plot_grid(plotlist = plot_list, ncol = 3)  # Adjust ncol for the number of columns you want

# Or using patchwork
library(patchwork)
combined_plot <- wrap_plots(plot_list, ncol = 3) +
  plot_annotation(title = "Plots of Estimated Proportions", 
                  theme = theme(plot.title = element_text(size = 14, face = "bold")))

print(combined_plot)
