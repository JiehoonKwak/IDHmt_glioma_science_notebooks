library(Seurat)
library(tidyverse)
library(here)
library(glue)
library(schard)

library(infercnv)
library(SCEVAN)
library(copykat)

# prepare gene annotations with gtf file with infercnvpy

# Prepare variables
query  <- "Mouse1"
n_jobs  <- 48
organism <- 'mouse'
genome <- ifelse(organism == 'mouse', 'mm10', 'hg20')

data_path <- here('notebook/tmp')
# infercnv reference cat
ref_group_names <- c("Microglia")

# Load data
obj  <- schard::h5ad2seurat(here(data_path, glue(query, '.h5ad')))
count <- GetAssayData(obj, layer = 'counts')
table(obj[[]]$cell_type)

obj@meta.data["cell_type"] |> rownames_to_column('index')  |> write_tsv(here(data_path, glue(query, '_cell_type.txt')), col_names = FALSE)

obj[["RNA"]]@meta.features |> rownames_to_column('index') %>% filter(complete.cases(.)) |> select(index, chromosome, start, end) |> write_tsv(here(data_path, glue(query, '_gene_annotation.txt')), col_names = FALSE)

# SCEVAN : output stored in './output' -> move to './scevan/{query}'
scevan_results <- pipelineCNA(count, sample = query, organism = organism, par_cores = n_jobs)
# scevan_results |> rownames_to_column('index') |> write_csv(here(data_path, glue(query, '_scevan.csv')) )

# CopyKat : output stored in './' -> move to './copykat/{query}'
copykat_results <- copykat(rawmat = count, KS.cut=0.1, sam.name = query, plot.genes="FALSE", genome = genome ,n.cores = n_jobs)
copykat_results$prediction |> write_csv(here(data_path, glue(query, '_copykat.csv')))

# InferCNV : output stored in ''
cnv <- CreateInfercnvObject(count, here(data_path, glue(query, '_gene_annotation.txt')), here(data_path, glue(query, '_cell_type.txt')),
                            ref_group_names=ref_group_names)

dir.create(here(data_path, glue('infercnv/', query)), recursive = TRUE, showWarnings = FALSE)

infercnv_results <- infercnv::run(cnv, cutoff=0.1, out_dir=here(data_path, glue('infercnv/', query)), denoise=TRUE, HMM=TRUE, num_threads = n_jobs)

# extract predictions
scevan_results # index + class
copykat_results$prediction # cell.names, copykat.pred