### SCANVI script
## conda install scvi-tools -c conda-forge
## conda install -c conda-forge -c bioconda scrublet
## seurat - anndata transformation:
## conda install -c conda-forge r-base r-seurat anndata scanpy
## pip install anndata2ri
## https://docs.scvi-tools.org/en/1.4.1/tutorials/notebooks/multimodal/scarches_scvi_tools.html
# use a reference and then work on query data

import sys
sys.stderr=open(snakemake.log[0], "a+")
sys.stdout=open(snakemake.log[0], "a+")

print("CONFIGURATION STEP")
# A. Parameters: 
# 1. Load libraries.
import os
import tempfile

import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scrublet as scr
import scvi
import seaborn as sns
import torch
import random

import scipy.sparse
import rpy2.robjects as robjects
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.conversion import localconverter

print("Last run with scvi-tools version:", scvi.__version__)

sc.set_figure_params(figsize=(6, 6), frameon=False)
sns.set_theme()
torch.set_float32_matmul_precision("high")

# %config InlineBackend.print_figure_kwargs={"facecolor": "w"}
# %config InlineBackend.figure_format="retina"

# 2. Folder configuration. 
ann_obj = snakemake.input["seurat_obj"] # 6_annotation rds object
outdir = snakemake.params["output_dir"]
random_seed = snakemake.params["random_seed"]
random.seed(random_seed)
np.random.seed(random_seed)
sc.settings.seed = random_seed
scvi.settings.seed = random_seed

folders = ["1_preprocessing", "2_normalization", "3_clustering", "4_annotation", "5_scRepertoire", "6_degs", "7_gs", "8_traj_in", "9_func_analysis", "10_RNAvelocity", "11_CellRank", "scANVI"]

input_dir =  os.path.join(outdir, folders[3])
output_dir = os.path.join(outdir, folders[3], folders[11])
os.makedirs(output_dir, exist_ok=True)

print("2. Folder paths were set.")
n_threads = snakemake.threads
mem_mb = snakemake.resources["mem_mb"]
walltime = snakemake.resources["walltime"]

sc.settings.n_jobs = n_threads
# Default scanpy save:
sc.settings.figdir = output_dir
sc.set_figure_params(frameon=False, dpi=150, facecolor='white')

def save_graph(file_name):
  path = os.path.join(output_dir, file_name)
  plt.savefig(path, bbox_inches='tight')
  plt.close()
  print(f"Saved at: {path}")
    
r('library(Seurat)')
r('library(Matrix)')

# Function to extract seurat original data slots for adata:
def get_sparse_matrix_from_r(r_obj_name, slot="counts", assay=None):
    if assay:
        cmd = f'GetAssayData({r_obj_name}, assay="{assay}", slot="{slot}")'
    else:
        cmd = f'GetAssayData({r_obj_name}, layer="{slot}")'
    
    r(f'temp_mat <- as({cmd}, "dgCMatrix")')
    
    data = np.array(r('temp_mat@x'))
    indices = np.array(r('temp_mat@i'))
    indptr = np.array(r('temp_mat@p'))
    shape = tuple(r('dim(temp_mat)')) # (Genes, Cells)
    
    # Reconstruct matrix for scipy
    mat_csc = scipy.sparse.csc_matrix((data, indices, indptr), shape=shape)
   
    return mat_csc.T

print("Configuration finished.")
print("\n")

print("PROCESSING STEP")
robjects.r.assign("ann_obj_r", ann_obj)
r_script = """
library(Seurat)
seu <- readRDS(ann_obj_r)
seu
"""
seurat_obj = r(r_script)

print("1. Seurat object was loaded.")

# Import expression matrix
X_sparse = get_sparse_matrix_from_r('seu', slot='counts')
cells = list(r('colnames(seu)'))
genes = list(r('rownames(seu)'))

# Import metadata
with localconverter(robjects.default_converter + pandas2ri.converter):
  meta_data = r('seu@meta.data')
  if not isinstance(meta_data, pd.DataFrame):
    meta_data = pandas2ri.rpy2py(meta_data)

# create adata
adata = sc.AnnData(X=X_sparse)
adata.obs_names = cells
adata.var_names = genes
adata.obs = meta_data

# warning: metadata objects might be a problem - clean
for col in adata.obs.columns:
    if adata.obs[col].dtype == 'object':
        try:
            adata.obs[col] = adata.obs[col].astype(str)
        except Exception:
            adata.obs[col] = adata.obs[col].apply(lambda x: str(x))
            
# Import UMAP & PCA
try:
  pca_coords = np.array(r('Embeddings(seu, reduction = "pca")'))
  umap_coords = np.array(r('Embeddings(seu, reduction = "umap")'))
  adata.obsm['X_pca'] = pca_coords
  adata.obsm['X_umap'] = umap_coords
  print("-- PCA & UMAP loaded.")
except Exception as e:
  print("Error when loading UMAP coordinates.")


adata.obs['seurat_clusters'] = adata.obs['seurat_clusters'].astype('category')

r_var_genes = r('VariableFeatures(object = seu)')
var_genes_list = list(r_var_genes)
adata.var['highly_variable'] = adata.var_names.isin(var_genes_list)

sc.pl.umap(adata, color="seurat_clusters", palette="tab20", show=False)
adata.write(os.path.join(input_dir, "seurat_annotation.h5ad"))

print("2. Seurat object converted to AnnData and properly parsed for scANVI analysis...")
# https://docs.scvi-tools.org/en/1.4.1/tutorials/notebooks/scrna/seed_labeling.html

# fulfill scType annotation (correct Unclassified cells)
adata.obs["seed_labels"] = adata.obs["scType_annotation_non_restrict"].astype(str)

adata.obs.loc[adata.obs["seed_labels"] == "Unclassified", "seed_labels"] = "Unknown"

print("scANVI modeling...")
unique_labels = adata.obs["seed_labels"].unique()
real_labels = [l for l in unique_labels if l != "Unknown"]

print(f"Found {len(real_labels)} distinct cell types: {real_labels}")
if len(real_labels) < 2:
  print("WARNING: Not enough cell types found to train a classifier (< 2). Skipping scANVI.")
  adata.obs["scANVI_annotation"] = adata.obs["scType_annotation_non_restrict"].astype(str)
  sc.pl.embedding(adata, basis="X_umap", color="scANVI_annotation", title="scANVI structure")
  save_graph("scANVI_predictions-seurat_umap.pdf")
  df_predict = pd.DataFrame(
    data = adata.obs["scANVI_annotation"].astype(str).values, 
    index = adata.obs_names, 
    columns = ["scANVI_annotation"]
  )
  tsv_path = os.path.join(input_dir, "scanvi_predicted_ann.tsv")
  df_predict.to_csv(tsv_path, sep='\t', header=True, index=True)
  robjects.globalenv['tsv_path_r'] = tsv_path
  r_script_save = """
  print("Exporting metadata to seurat object...")
  preds <- read.table(tsv_path_r, sep="\t", header=TRUE, row.names=1)

  seu <- AddMetaData(object = seu, metadata = preds)

  # verify column creation
  print(table(seu@meta.data$scANVI_annotation))

  output_path <- paste0(dirname(ann_obj_r), "/seurat_annotated_scanvi.rds")
  saveRDS(seu, file = output_path)
  """

  robjects.r(r_script_save)
  sys.stderr.close()
  sys.stdout.close()
  sys.exit(0)
    
scvi.model.SCVI.setup_anndata(adata, batch_key=None, labels_key="seed_labels")
scvi_model = scvi.model.SCVI(adata, n_latent=30, n_layers=2)
scvi_model.train(100)

scanvi_model = scvi.model.SCANVI.from_scvi_model(scvi_model, "Unknown")
scanvi_model.train(25)

print("scANVI ppredicting...")
SCANVI_LATENT_KEY = "X_scANVI"
SCANVI_PREDICTIONS_KEY = "C_scANVI"

adata.obsm[SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation(adata)
adata.obs[SCANVI_PREDICTIONS_KEY] = scanvi_model.predict(adata)

adata.obsm["X_umap_seurat"] = adata.obsm["X_umap"].copy() # scANVI UMAP will replace seurat UMAP

sc.pp.neighbors(adata, use_rep=SCANVI_LATENT_KEY, key_added="neighbors_scanvi")
sc.tl.umap(adata, neighbors_key="neighbors_scanvi")
adata.obsm["X_umap_scanvi"] = adata.obsm["X_umap"].copy()

sc.pl.embedding(adata, basis="X_umap_seurat", color=SCANVI_PREDICTIONS_KEY, title="Estructura Seurat")
save_graph("scANVI_predictions-seurat_umap.pdf")

sc.pl.embedding(adata, basis="X_umap_scanvi", color=SCANVI_PREDICTIONS_KEY, title="Estructura scANVI")
save_graph("scANVI_predictions-scANVI_umap.pdf")

print("3- scanVI analysis completed")

adata.write(os.path.join(output_dir, "seurat_scanvi.h5ad"))
print("4- AnnData object saved")

df_predict = pd.DataFrame(
    data = adata.obs[SCANVI_PREDICTIONS_KEY].astype(str).values, 
    index = adata.obs_names, 
    columns = ["scANVI_annotation"]
)
tsv_path = os.path.join(input_dir, "scanvi_predicted_ann.tsv")
df_predict.to_csv(tsv_path, sep='\t', header=True, index=True)

robjects.globalenv['tsv_path_r'] = tsv_path

r_script_save = """
print("Exporting metadata to seurat object...")

preds <- read.table(tsv_path_r, sep="\t", header=TRUE, row.names=1)

seu <- AddMetaData(object = seu, metadata = preds)

# verify column creation
print(table(seu@meta.data$scANVI_annotation))

output_path <- paste0(dirname(ann_obj_r), "/seurat_annotated_scanvi.rds")
saveRDS(seu, file = output_path)

"""

robjects.r(r_script_save)

print("5- Annotation exported to seurat_annotated.rds")

sys.stderr.close()
sys.stdout.close()

