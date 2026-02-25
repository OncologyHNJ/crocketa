# conda install -c conda-forge cellrank
# conda install -c conda-forge r-base r-seurat anndata scanpy
# pip install anndata2ri

import sys
sys.stderr=open(snakemake.log[0], "a+")
sys.stdout=open(snakemake.log[0], "a+")
print("CONFIGURATION STEP")
# A. Parameters: 
# 1. Load libraries.
import os
import random
import numpy as np
import matplotlib
matplotlib.use("Agg") # avoid internal errors related to graphical devices
import cellrank as cr
import anndata2ri
import scanpy as sc
import pandas as pd
import scipy.sparse
import rpy2.robjects as robjects
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.conversion import localconverter
import scvelo as scv
scv.settings.verbosity = 3
scv.settings.set_figure_params("scvelo")
cr.settings.verbosity = 2
import warnings
import matplotlib.pyplot as plt
from cellrank.estimators import GPCCA

# 2. Folder configuration. 
vel_obj = snakemake.input["seurat_obj"]
outdir = snakemake.params["output_dir"]
random_seed = snakemake.params["random_seed"]
random.seed(random_seed)
np.random.seed(random_seed)
sc.settings.seed = random_seed
folders = ["1_preprocessing", "2_normalization", "3_clustering", "4_annotation", "5_scRepertoire", "6_degs", "7_gs", "8_traj_in", "9_func_analysis", "10_RNAvelocity", "11_CellRank"]
input_dir =  os.path.join(outdir, folders[9])
output_dir = os.path.join(outdir, folders[10])
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
robjects.r.assign("vel_obj_r", vel_obj)
r_script = """
library(Seurat)
seu <- readRDS(vel_obj_r)
vel_data <- Tool(object = seu, slot = "RunVelocity")
seu[['velocity_results']] <- CreateAssayObject(counts = vel_data$deltaE)
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

# Import UMAP & PCA
try:
  pca_coords = np.array(r('Embeddings(seu, reduction = "pca")'))
  umap_coords = np.array(r('Embeddings(seu, reduction = "umap")'))
  adata.obsm['X_pca'] = pca_coords
  adata.obsm['X_umap'] = umap_coords
  print("-- PCA & UMAP loaded.")
except Exception as e:
  print("Error when loading UMAP coordinates.")
  
# Import velocity
try:
  vel_matrix = get_sparse_matrix_from_r('seu', slot='counts', assay='velocity_results')
  vel_genes = list(r('rownames(seu[["velocity_results"]])'))
  df_vel = pd.DataFrame.sparse.from_spmatrix(
    vel_matrix, 
    columns=vel_genes, 
    index=adata.obs_names
  )
  df_vel_aligned = df_vel.reindex(columns=adata.var_names, fill_value=0)
  r_var_genes = r('VariableFeatures(object = seu)')
  var_genes_list = list(r_var_genes)
  adata.var['highly_variable'] = adata.var_names.isin(var_genes_list)
  adata.layers['velocity'] = scipy.sparse.csr_matrix(df_vel_aligned.values)
  adata.var['velocity_genes'] = adata.var_names.isin(vel_genes)
  print("-- Velocity loaded.")
  
except Exception as e:
  print(f"Error when loading Velocity: {e}")

# --- Verify ---
print("adata dims:", adata.shape)
print("velocity dims:", adata.layers['velocity'].shape)

# Import spliced/unspliced layers
def import_and_align_layer(r_assay_name, layer_name_in_scanpy):
	print(f"-- Processing layer: {layer_name_in_scanpy} ---")
	raw_matrix = get_sparse_matrix_from_r('seu', assay=r_assay_name, slot='counts')
	raw_genes = list(r(f'rownames(GetAssayData(seu, assay="{r_assay_name}"))'))
	print(f"-- Original matrix loaded: {raw_matrix.shape}")
	df_temp = pd.DataFrame.sparse.from_spmatrix(
		raw_matrix, 
		columns=raw_genes, 
		index=adata.obs_names
	)
	df_aligned = df_temp.reindex(columns=adata.var_names, fill_value=0)
	adata.layers[layer_name_in_scanpy] = scipy.sparse.csr_matrix(df_aligned.values)
	print(f"-- Succeeded: Layer '{layer_name_in_scanpy}' aligned to {df_aligned.shape}")

# a- Import SPLICED
try:
  import_and_align_layer('spliced', 'spliced')
except Exception as e:
  print(f"Error when importing spliced layer: {e}")

# b- Import UNSPLICED
try:
  import_and_align_layer('unspliced', 'unspliced')
except Exception as e:
  print(f"Error when importing unspliced layer: {e}")

print("\nData layers contained in adata: ", adata.layers.keys())

adata.obs['seurat_clusters'] = adata.obs['seurat_clusters'].astype('category')

sc.pl.umap(adata, color="seurat_clusters", palette="tab20", show=False)

adata.write(os.path.join(input_dir, "seurat_velocity.h5ad"))

# keep only relevant genes: velocity genes + highly variable genes
adata.var['velocity_genes'] = np.array(np.abs(adata.layers['velocity']).sum(0)).flatten() > 0

genes_to_keep = adata.var['velocity_genes'] | adata.var['highly_variable']
adata = adata[:, genes_to_keep].copy()

def make_dense_and_clean(matrix):
  # if sparsed matrix, convert to dense
  if hasattr(matrix, "toarray"):
    matrix = matrix.toarray()
  # match object type (float32)
  matrix = np.array(matrix, dtype=np.float32)
  # clean NAs/Inf
  matrix = np.nan_to_num(matrix, nan=0.0, posinf=0.0, neginf=0.0)
  return matrix

adata.X = make_dense_and_clean(adata.X)
adata.layers['velocity'] = make_dense_and_clean(adata.layers['velocity'])
adata.layers['spliced'] = make_dense_and_clean(adata.layers['spliced'])
adata.layers['unspliced'] = make_dense_and_clean(adata.layers['unspliced'])

adata.write(os.path.join(input_dir, "seurat_velocity_filtered.h5ad"))

print("2. Seurat object converted to AnnData and properly parsed for analysis...")

print("-- Recalculate nearest neighbors")
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30, use_rep='X_pca', random_state=random_seed)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.pl.proportions(adata, show = False)
save_graph("1-Proportions_splicing.png")
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, basis='umap', color='seurat_clusters', legend_loc='right margin', show=False)
save_graph("2-Velocity_Stream.png")

vk = cr.kernels.VelocityKernel(adata)
vk.compute_transition_matrix()

# combine velocity with connectivity kernel (gene expression similarity)
ck = cr.kernels.ConnectivityKernel(adata)
ck.compute_transition_matrix()

combined_kernel = 0.6 * vk + 0.4 * ck # ensure general analysis: 60% velocity, 40% gene expression

vk.plot_projection(basis="umap", color="seurat_clusters", legend_loc='right margin', show = False)
save_graph("3-Projection_kernel.png")

## Infer starting/ending clusters mathematically:
g = GPCCA(combined_kernel)
g.compute_schur(n_components=20)
g.plot_spectrum(real_only = True, legend_loc='upper right' )
save_graph("4.a-Eigengap_Spectrum.png")
try:
  print("-- Automatically infer Macrostates...")
  g.compute_macrostates(
    cluster_key='seurat_clusters'
  ) 
except Exception as e:
  print(f"Error extracting macrostates: {e}")
  print("Try again with default n_states=3...")
  g.compute_macrostates(cluster_key='seurat_clusters', n_states=3)

g.plot_macrostates(which='all', basis='umap', title='Detected macrostates', legend_loc='right margin', show=False)
save_graph("4.b-Macrostates_Map.png")

g.predict_terminal_states(method = "eigengap")
g.predict_initial_states(allow_overlap = True)
g.compute_fate_probabilities()

g.plot_macrostates(which='initial', basis='umap', discrete = True, title='Infered terminal macrostates', legend_loc='right margin', show = False)
save_graph("5-Terminal_States.png")

g.plot_macrostates(which='terminal', basis='umap', discrete = True, title='Infered initial macrostates', legend_loc='right margin', show = False)
save_graph("6-Initial_States.png")

g.plot_fate_probabilities(basis='umap', legend_loc='right margin',show=False)
save_graph("7-Fate_Probabilities.png")

print("Extract random walk:")
try:
  initial_clusters = g.initial_states.cat.categories.tolist()
  print(f"-- Initial states detected: {initial_clusters}")
  start_config = {'seurat_clusters': initial_clusters}
  vk.plot_random_walks(
    start_ixs=start_config,
    max_iter=200,
    seed=random_seed,
    basis='umap',
    title=f"Trajectories from {initial_clusters}",
    legend_loc='right margin'
  )
  save_graph("8-Random_Walks_Trajectory.png")
except AttributeError:
  print("Initial states from object could not be infered.")
except Exception as e:
  print(f"Error when drawing: {e}")


print("3- CellRank analysis completed: CellRank + Velocity")

vk.write_to_adata()
adata.write(
    os.path.join(output_dir, "seurat_cellrank.h5ad")
)


print("4- AnnData object saved")

sys.stderr.close()
sys.stdout.close()
