devtools::install_github('satijalab/seurat-data')
library(SeuratData)
library(Seurat)

InstallData("panc8")
data("panc8")

# Enforce Seurat object version compatibility, if needed
panc8 <- UpdateSeuratObject(panc8)

panc8 <- panc8 |> NormalizeData() |>
  FindVariableFeatures() |>
  ScaleData() |> RunPCA(npcs=20)



# We can calculate metrics for batch mixing and cell type separation on this 
# unintegrated object, by specifying the metadata columns that contain batch 
# and cell type information ("tech" and "celltype" respectively):

library(scIntegrationMetrics)
metrics <- getIntegrationMetrics(panc8, meta.label = "celltype",
                                 meta.batch = "tech",
                                 iLISI_perplexity = 20)
unlist(metrics)


# We can apply an integration method such as STACAS to mitigate batch effects:
remotes::install_github("carmonalab/STACAS")
library(STACAS)
panc8.list <- SplitObject(panc8, split.by = "tech")

panc8.stacas <- Run.STACAS(panc8.list)

# We can then calculate metrics after integration, and verify whether batch 
# mixing and cell type separation were improved upon integration:
  
metrics.integrated <- getIntegrationMetrics(panc8.stacas, 
                                            meta.label = "celltype",
                                            meta.batch = "tech",
                                            iLISI_perplexity = 20)
unlist(metrics.integrated)


#############################################################

# test integration metrics in merges samples but without joining layers
samples.unintegrated <- getIntegrationMetrics(samples, 
                                            meta.label = "unintegrated_clusters",
                                            meta.batch = "sample.id",
                                            iLISI_perplexity = 20)

unlist(samples.unintegrated)
# start 250520 1:34pm to error
# start 250520 2:03pm to 2:07pm, took 4 min


# test integration metrics in merges samples but without joining layers
samples.integrated <- getIntegrationMetrics(sampleswi, 
                                              meta.label = "seurat_clusters",
                                              meta.batch = "sample.id",
                                              iLISI_perplexity = 20)

unlist(samples.integrated)
# start 250520 3:09pm




samples.unintegrated <- getIntegrationMetrics(samples, 
                                              meta.label = "unintegrated_clusters",
                                              meta.batch = "sample.id",
                                              iLISI_perplexity = 20)

unlist(samples.unintegrated)

#################################################################################

samples.integrated <- getIntegrationMetrics(sampleswi, 
                                            meta.label = "seurat_clusters",
                                            meta.batch = "sample.id",
                                            iLISI_perplexity = 20)

unlist(samples.integrated)


#################################################################################

samples.harmony <- getIntegrationMetrics(samples_harmony, 
                                            meta.label = "harmony.cluster",
                                            meta.batch = "sample.id",
                                            iLISI_perplexity = 20)

unlist(samples.harmony)











































