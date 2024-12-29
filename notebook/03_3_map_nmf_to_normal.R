library(Seurat)
library(tidyverse)
library(here)
library(schard)
options(future.globals.maxSize = 100 * 1024^3) 

tumor  <-  h5ad2seurat(here('output/9_annotated_subcluster_embedding_tumor_nmf_added.h5ad'))
tumor <- SCTransform(tumor, vars.to.regress = "sample_id") |> 
    RunPCA(verbose = FALSE) |> 
    RunUMAP(dims = 1:30, verbose = FALSE)
DimPlot(tumor, group.by = c("cell_type", "sample_id"))


adata_f  <-  h5ad2seurat(here('public/allen_fetal_nmf_filtered.h5ad'))
adata_f <- SCTransform(adata_f, vars.to.regress = "sample", verbose = FALSE) |> 
    RunPCA(verbose = FALSE) |> 
    RunUMAP(dims = 1:30, verbose = FALSE)
DimPlot(adata_f, reduction = "umap", group.by = c("cell_type", "sample"), label = TRUE)


adata_a  <-  h5ad2seurat(here('public/allen_fetal_nmf_filtered.h5ad'))
adata_a
adata_a <- SCTransform(adata_a, vars.to.regress = "sample", verbose = FALSE) |> 
    RunPCA(verbose = FALSE) |> 
    RunUMAP(dims = 1:30, verbose = FALSE)
DimPlot(adata_a, reduction = "umap", group.by = c("cell_type", "sample"), label = TRUE)