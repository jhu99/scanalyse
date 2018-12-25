import scanpy.api as sc
import numpy as np
adata = sc.read("./data/scRNAseq_CountTable.csv", first_column_names=True)
adata = adata.transpose()

print(adata)
sc.pp.filter_genes(adata, min_counts=1)
print(adata)
sc.pp.filter_cells(adata, min_counts=1)
print(adata)

adata.raw = adata.copy()

sc.pp.normalize_per_cell(adata)
print(adata.obs.n_counts)
adata.obs['size_factors'] = adata.obs.n_counts / np.median(adata.obs.n_counts)
print(sum(adata.X))
sc.pp.log1p(adata)
print(sum(adata.X))
sc.pp.scale(adata)
print(sum(adata.X))
