import sys
import os
import scanpy as sc
import scvelo as scv
import pandas as pd
import numpy as np
import loompy

scv.settings.set_figure_params('scvelo')

# merge loom files
ctrl_fn = "/Users/kriemo/Projects/sc_repos/vanotterloo/data/cellranger/1_control/velocyto/1_control.loom"
mut_fn = "/Users/kriemo/Projects/sc_repos/vanotterloo/data/cellranger/2_mutant/velocyto/2_mutant.loom"


loom_dir = "objects"
if not os.path.exists(loom_dir):
    os.makedirs(loom_dir)

if not os.path.isfile(os.path.join(loom_dir, "combined.loom")):
  loompy.combine([ctrl_fn, mut_fn], os.path.join(loom_dir, "combined.loom"))
 
 
# load
adata = scv.read(os.path.join(loom_dir, "combined.loom"), cache=True)

adata.var_names_make_unique()

mdata = pd.read_csv("results/tables/ameloblast_metadata_2019_07_07.tsv.gz", sep="\t")

# get cell ids to match loom object
new_ids = []
for cell in mdata["cell"]:
  fields = cell.split("_")
  sample_id = fields[0] + "_" + fields[1]
  bc = fields[2] + "x"
  new_id = sample_id + ":" + bc
  new_ids.append(new_id)
  
mdata = mdata.assign(new_id = new_ids)

mdata = mdata[mdata.new_id.isin(list(adata.obs.index.values))]

# reorder cell ids to match loom object
cids = pd.DataFrame({'new_id' : adata.obs.index.values})
mdata = pd.merge(cids, mdata, how = 'left', on = 'new_id')
mdata = mdata.dropna()

#only keep cells found in seurat data
keep_idx = [x in list(mdata["new_id"]) for x in list(adata.obs.index.values)]

adata = adata[keep_idx, :]
adata = adata.copy()

#add cluster annotations

adata.obs['clusters'] = np.array([str(x) for x in mdata["clusters"]])
adata.obs['sample'] = np.array(mdata["sample"])

#add tSNE projections
umap_mat = np.column_stack((np.array(mdata["UMAP_1"]), np.array(mdata["UMAP_2"])))
adata.obsm['X_umap'] = umap_mat

scv.utils.show_proportions(adata)

scv.pp.filter_genes(adata, min_counts=20, min_counts_u=10)
scv.pp.normalize_per_cell(adata)
scv.pp.filter_genes_dispersion(adata, n_top_genes=3000)
scv.pp.log1p(adata)


scv.pp.moments(adata, n_pcs=20, n_neighbors=30)
scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata)
scv.tl.velocity_embedding(adata, basis='umap')
