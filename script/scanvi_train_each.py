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
# 1. setup anndata
adata = sc.read_h5ad('../output/7_annotated_cnv.h5ad')
sc.pp.highly_variable_genes(adata, n_top_genes = 3000, subset = True, layer = 'counts', flavor = 'seurat_v3', batch_key = batch_key) 
model = scvi.model.SCVI.load('../output/model_hvg_all.model', adata)
print('scVI model loaded')

# 2. train scANVI
scanvi_model = scvi.model.SCANVI.from_scvi_model(scvi_model= model,adata = adata, labels_key = 'each_cell_type', unlabeled_category='unlabelled')

scanvi_model.train(max_epochs = 30, early_stopping = True)
print('scANVI model trained')

with mplscience.style_context():
    y = scanvi_model.history['reconstruction_loss_validation']['reconstruction_loss_validation'].min()
    plt.plot(scanvi_model.history['reconstruction_loss_validation']['reconstruction_loss_validation'], label = 'validation')
    plt.plot(scanvi_model.history['reconstruction_loss_train']['reconstruction_loss_train'], label = 'train')

    plt.axhline(y, c = 'k')
    plt.legend()
    plt.savefig('../plot/scanvi_loss.png')
    plt.close()
    
scanvi_model.save('../output/model_scanvi')
print('scANVI model saved')

# 4. save results
adata = sc.read_h5ad('../output/7_annotated_cnv.h5ad')
adata.obsm["X_scANVI"] = scanvi_model.get_latent_representation()
sc.pp.neighbors(adata, use_rep = 'X_scANVI')
sc.tl.umap(adata)

adata.write_h5ad('../output/8_scANVI_each.h5ad')
print('scANVI results saved')