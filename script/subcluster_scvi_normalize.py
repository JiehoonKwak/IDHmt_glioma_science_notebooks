import os
import warnings
warnings.filterwarnings("ignore")
import scanpy as sc
import numpy as np
import scvi
import torch
import matplotlib.pyplot as plt
import mplscience

np.random.seed(777)
torch.manual_seed(777)
scvi.settings.seed = 777
batch_key = 'sample_id'

adata = sc.read_h5ad('../output/annotated_subcluster.h5ad')
scvi.model.SCVI.setup_anndata(adata, layer = 'counts', batch_key = batch_key, continuous_covariate_keys=['pct_counts_mt', 'pct_counts_ribo'])
model_scvi = scvi.model.SCVI(adata, n_layers=2)
model_scvi.train(early_stopping=True)
model_scvi.save('../output/model_subcluster_all')
adata.X = adata.layers['counts'].toarray()
adata.obsm['X_scVI'] = model_scvi.get_latent_representation()
adata.layers['scvi_normalized'] = model_scvi.get_normalized_expression(library_size=1e6) # tpm

adata.write('../output/subcluster.h5ad')

adata.X = adata.layers['scvi_normalized'].copy()
adata.write_h5ad('../output/subcluster_scvi.h5ad')

sc.pp.highly_variable_genes(adata, n_top_genes = 2000, subset = True, layer = 'counts', flavor = 'seurat_v3', batch_key = 'sample_id') 