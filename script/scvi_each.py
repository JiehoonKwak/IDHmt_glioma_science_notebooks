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

# 1. setup anndata
batch_key = 'sample_id'
adata = sc.read_h5ad('../output/7_annotated_cnv.h5ad')
sc.pp.highly_variable_genes(adata, n_top_genes = 3000, subset = True, layer = 'counts', flavor = 'seurat_v3', batch_key = batch_key) 

# 2. train scVI
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=batch_key, continuous_covariate_keys=['pct_counts_mt', 'pct_counts_ribo'])

model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="zinb", dropout_rate = 0.5)

model.train(early_stopping=True, datasplitter_kwargs={"drop_last": True}, plan_kwargs={"lr": 0.0027},)

with mplscience.style_context():
    y = model.history['reconstruction_loss_validation']['reconstruction_loss_validation'].min()
    plt.plot(model.history['reconstruction_loss_validation']['reconstruction_loss_validation'], label = 'validation')
    plt.plot(model.history['reconstruction_loss_train']['reconstruction_loss_train'], label = 'train')

    plt.axhline(y, c = 'k')
    plt.legend()
    plt.savefig('../plot/scvi_loss.png')
    plt.close()
    
model.save('../output/model_cnv')

# 3. train scANVI
scanvi_model = scvi.model.SCANVI.from_scvi_model(
    model,
    adata=adata,
    labels_key="each_cell_type",
    unlabeled_category="Unknown",
)
scanvi_model.train(max_epochs=20, early_stopping=True)
with mplscience.style_context():
    y = scanvi_model.history['reconstruction_loss_validation']['reconstruction_loss_validation'].min()
    plt.plot(scanvi_model.history['reconstruction_loss_validation']['reconstruction_loss_validation'], label = 'validation')
    plt.plot(scanvi_model.history['reconstruction_loss_train']['reconstruction_loss_train'], label = 'train')

    plt.axhline(y, c = 'k')
    plt.legend()
    plt.savefig('../plot/scanvi_loss.png')
    plt.close()
    
scanvi_model.save('../output/model_scanvi')
    
# 4. save results
adata = sc.read_h5ad('../output/7_annotated_cnv.h5ad')
adata.obsm['X_scVI_re'] = model.get_latent_representation()
adata.obsm['X_scANVI'] = scanvi_model.get_latent_representation()

adata.write_h5ad('../output/8_scANVI_each.h5ad')
