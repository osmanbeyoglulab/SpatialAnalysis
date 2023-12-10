library(Seurat)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)

setwd("~/Desktop/Pitt/OV")
source('aux_spatial.R')

n = 11
file.dirs = paste0('seurat_raw/GSE189843_Stur/II214', seq(72,83))
file.dirs[8:12] = file.dirs[c(9:12,8)]
sample.ids = c(paste0('ER_', 1:6), paste0('PR_', 1:5))
group.ids = c(rep('ER', 6), rep('PR', 5))

# ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====

# Loading
ov = list()
for (i in seq(n)){
  ov[[sample.ids[i]]] = Load10X_Spatial_v2(file.dirs[i], slice = sample.ids[i])
  ov[[sample.ids[i]]]$Sample = sample.ids[i]
  ov[[sample.ids[i]]]$Group = group.ids[i]
}

# Normalization
for (sample in sample.ids){
  ov[[sample]] = SCTransform(ov[[sample]], assay = "Spatial") %>%
    RunPCA(assay = "SCT", verbose = T) %>%
    FindNeighbors(reduction = "pca", dims = 1:30) %>%
    FindClusters(resolution = 2, verbose = T) %>%
    RunUMAP(reduction = "pca", dims = 1:30)
}

# Spatially variable features
for (sample in sample.ids){
  ov[[sample]] = FindSpatiallyVariableFeatures(
    ov[[sample]], assay = "SCT", 
    features = VariableFeatures(ov[[sample]])[1:1000],
    selection.method = c("markvariogram", "moransi"),
    verbose = T)
}

saveRDS(ov, sprintf("seurat_spatial/gse189843_stur_sct%s.rds", n))
