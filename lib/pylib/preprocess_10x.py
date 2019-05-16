import scipy as sci
import numpy as np
from scipy.sparse import csr_matrix, csc_matrix
import scanpy as sc
import anndata
from anndata import h5py, logging
from anndata import AnnData
import pandas as pd
import matplotlib
import matplotlib.pyplot as pl
import seaborn as sns

matplotlib.use('Agg')
sns.set(style='white', rc={'figure.figsize':(8,6), 'figure.dpi':150})

def getAnnData(input_file):
	h5 = h5py.File(input_file,'r')
	data = h5['/GRCh38/data']
	indices = h5['/GRCh38/indices']
	barcodes = h5['/GRCh38/barcodes']
	indptr = h5['/GRCh38/indptr']
	genes = h5['/GRCh38/genes']
	gene_names = h5['/GRCh38/gene_names']
	shape = h5['/GRCh38/shape']

	X = csr_matrix((data,indices,indptr),shape=(shape.value[1],shape.value[0]),dtype='float32')
	adata = AnnData(X,
		obs=pd.DataFrame(index=barcodes.value),
		var=pd.DataFrame(index=genes.value),
		dtype=X.dtype.name,
		filemode=True)
	return adata

def getNormAnnData(input_file,normalize_type):
		h5 = h5py.File(input_file,'r')
		data = h5['/GRCh38/data']
		indices = h5['/GRCh38/indices']
		barcodes = h5['/GRCh38/barcodes']
		indptr = h5['/GRCh38/indptr']
		genes = h5['/GRCh38/genes']
		gene_names = h5['/GRCh38/gene_names']
		shape = h5['/GRCh38/shape']

		if normalize_type == "rank":
				zero_value = h5['/GRCh38/zero_value']
		else:
				size_factor = h5['/GRCh38/size_factor'][:]

		X = csr_matrix((data,indices,indptr),shape=(shape.value[1],shape.value[0]),dtype='float32')

		if	normalize_type == "rank":
				X = X.toarray()
				for i in range(shape.value[1]):
						for j in range(shape.value[0]):
								if abs(X[i][j] - 0) < 1e-5:
										X[i][j] = zero_value[i]

		adata = AnnData(X,
			 obs=pd.DataFrame(index=barcodes.value),
			 var=pd.DataFrame(index=genes.value),
			 dtype=X.dtype.name,
			 filemode=True)

		if normalize_type == "rank":
				return adata
		else:
				return adata,size_factor
 
def getAnnData_10x_h5(input_file):
	adata = sc.read_10x_h5(input_file,"GRCh38")
	return adata

def getAnnData_10x_mtx(input_file):
	adata = sc.read_10x_mtx(input_file)
	return adata

## reorder the variables in a desired list
def pre_process_input_data(adata,gene_file,format_type="10x_mtx", mode="test"):
	rownames = adata.obs_names.values
	colnames = adata.var['gene_ids'].values
	X2=pd.DataFrame(adata.X.todense(), index=rownames, columns=colnames)
	
	gene_input = pd.read_csv(gene_file,index_col=1)
	
	common_ind = pd.Index(colnames).intersection(gene_input.index.values)
	left_ind = gene_input.index.difference(common_ind)
	X3=X2[common_ind]
	for ind in left_ind:
		X3.insert(0,ind,0.01)
	X3=X3.reindex(columns=gene_input.index.values)
	X=sci.sparse.csr_matrix(X3,dtype='float32')
	adata = AnnData(X,
				obs=pd.DataFrame(index=rownames),
				var=pd.DataFrame(index=gene_input.index.values),
				dtype=X.dtype.name,
				filemode=True)
	if mode=="test":
		sc.pp.normalize_per_cell(adata)
		# used for loss evaluation
		adata.obs['size_factors'] = adata.obs.n_counts/ 10e4
		# np.median(adata.obs.n_counts)
		sc.pp.log1p(adata)
		sc.pp.scale(adata, zero_center=False)
		
	return adata

def filter_genes(adata,ntg=None,min_counts=None,min_percentage=None,flavor=None,method="HVG"):
	if method=='HVG':
		sc.pp.highly_variable_genes(adata,n_top_genes=ntg,flavor=flavor,inplace=True)
	elif method=='minexp':
		from math import ceil
		if min_counts is not None:
			sc.pp.filter_genes(adata,min_counts=min_counts)
		if min_percentage is not None:
			sc.pp.filter_genes(adata,min_cells=ceil(min_percentage*adata.shape[0]))
	elif method=='random':
		select_random()
	else:
		raise ValueError('`method` needs to be \'HVG\' or \'random\' or \'dropout\'')
	return adata
	
def filter_cells(adata,min_counts=None,min_percentage=None):
	from math import ceil
	if min_counts is not None:
		sc.pp.filter_cells(adata,min_counts=min_counts)
	if min_percentage is not None:
		sc.pp.filter_cells(adata,min_genes=ceil(min_percentage*adata.shape[1]))

def select_random():
	raise ValueError("wait for update")
	
def read_10x_data(input_file,format_type='10x_h5',backed=None):
	if format_type=='10x_h5':
		adata = sc.read_10x_h5(input_file)
	elif format_type=='10x_mtx':
		adata = sc.read_10x_mtx(input_file)
	elif format_type=='10x_h5ad':
		adata = sc.read_h5ad(input_file,backed=backed)
	elif format_type=="10x_csv":
		adata = sc.read_csv(input_file)
	else:
		raise ValueError('`format` needs to be \'10x_h5\' or \'10x_mtx\'')
	
	adata.var_names_make_unique()
	adata.obs_names_make_unique()
	return adata

def filter_basic(adata):
	sc.pp.filter_genes(adata,min_counts=500)
	sc.pp.filter_cells(adata,min_genes=200)
	mito_genes = adata.var_names.str.startswith('MT-')
	adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
	adata.obs['n_counts'] = adata.X.sum(axis=1).A1
	return adata
	
def add_annotation(adata,f):
	ann = pd.read_csv(f,sep='\t',header=0)
	adata.obs = adata.obs[['louvain','leiden']]
	# adata.obs = adata.obs.join(ann,lsuffix='_l',rsuffix='_r')
	
def recipe_zheng(adata,n_top_genes=1000, filter_disp_genes=True, path=None):
	if filter_disp_genes:
		sc.pp.filter_genes(adata,min_counts=1)
		sc.pp.filter_cells(adata,min_genes=200)
		mito_genes = adata.var_names.str.startswith('MT-')
		adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
		adata.obs['n_counts'] = adata.X.sum(axis=1).A1
		adata = adata[adata.obs['n_genes'] < 2500, :]
		adata = adata[adata.obs['percent_mito'] < 0.05, :]
		sc.pp.normalize_per_cell(adata,key_n_counts='n_counts_all',counts_per_cell_after=1e6)
		filter_result = sc.pp.filter_genes_dispersion(adata.X, flavor='cell_ranger', n_top_genes=n_top_genes, log=False)
		gene_list = adata.var_names[filter_result.gene_subset]
		gene_list = gene_list.union(pd.read_csv("./data/marker_genes.csv").x)
	else:
		gene_list=pd.read_csv(path+"genelist.csv",index_col=0).index
	adata = adata[:,gene_list]
	sc.pp.normalize_per_cell(adata,counts_per_cell_after=1e6)
	# used for loss evaluation
	if filter_disp_genes:
		adata.obs['size_factors'] = (adata.obs.n_counts_all / 1e6)*(adata.obs.n_counts/1e6)
	else:
		adata.obs['size_factors'] = adata.obs.n_counts/1e6
	sc.pp.log1p(adata)
	# sc.pp.scale(adata)
	if filter_disp_genes:
		adata.var.to_csv(path+'genelist.csv')
	return adata

def recipe_seurat(adata, log=True):
	sc.pp.recipe_seurat(adata)
	
def write2mtx(adata, path):
	gn=adata.var
	bc=adata.obs
	gn.to_csv(path+"genes.tsv",sep="\t",index=True,header=False)
	bc.to_csv(path+"barcodes.tsv",sep="\t",index=True,header=False)
	sci.io.mmwrite(path+"matrix.mtx", adata.X.transpose().astype(int))
	
def normalization(adata,path,filter_disp_genes=True):
	adataraw=adata.copy()
	adata=recipe_zheng(adata,path=path,filter_disp_genes=filter_disp_genes)
	cell_idx=adata.obs_names
	gene_idx=adata.var_names
	adata_subset = adataraw[:,gene_idx]
	adata_subset = adata_subset[cell_idx,:]
	adata.raw=adata_subset.copy()
	return adata
	
def nn_embedding(adata,path):
	sc.pp.neighbors(adata)
	sc.tl.louvain(adata)
	sc.tl.paga(adata)
	sc.pl.paga(adata, plot=False)
	sc.tl.umap(adata,init_pos="paga",min_dist=0.1)
	sc.tl.tsne(adata, n_jobs=12)
	adata.write(path+"ica_clusters.h5ad",compression='gzip')
	
def plotEmbedding(adata,path,color_col="louvain",ncol=10):
	sc.pl.umap(adata,color=color_col,show=False)
	pl.title("")
	pl.legend(loc=3,fontsize=6,mode="expand",bbox_to_anchor=(0.0, 1.01, 1, 0.2),ncol=ncol)
	pl.savefig(path+color_col+"_umap.png")
	pl.close()
	sc.pl.tsne(adata,color=color_col,show=False)
	pl.title("")
	pl.legend(loc=3,fontsize=6,mode="expand",bbox_to_anchor=(0.0, 1.01, 1, 0.2),ncol=ncol)
	pl.savefig(path+color_col+"_tsne.png")
	pl.close()

def plotTest(adata,path,t=None):
	if t is None:
		for t in range(10):
			adatat=adata[adata.obs['test']==t,:]
			tpath=path+"plot"+str(t)
			plotEmbedding(adatat,path=tpath,color_col='test')
	else:
		adatat=adata[adata.obs['test']==t,:]
		tpath=path+"plot"+str(t)
		plotEmbedding(adatat,path=tpath,color_col='test')
		
