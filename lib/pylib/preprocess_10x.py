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
	
def read_10x_data(input_file,format_type='10x_h5',backed=None,transpose=False,sparse=False):
	if format_type=='10x_h5':
		adata = sc.read_10x_h5(input_file)
	elif format_type=='10x_mtx':
		adata = sc.read_10x_mtx(input_file)
	elif format_type=='10x_h5ad':
		adata = sc.read_h5ad(input_file,backed=backed)
	elif format_type=="10x_csv":
		adata = sc.read_csv(input_file)
	elif format_type=="10x_txt":
		adata = sc.read_csv(input_file,delimiter="\t")
	else:
		raise ValueError('`format` needs to be \'10x_h5\' or \'10x_mtx\'')
	
	if transpose:
		adata = adata.transpose()
	if sparse:
		adata.X = csr_matrix(adata.X,dtype='float32')
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
	
def recipe_zheng(adata,n_top_genes=1000, filter_disp_genes=True, path=None, hg=True, ncounts=1e6, scale=False):
	if filter_disp_genes:
		sc.pp.filter_genes(adata,min_counts=1)
		sc.pp.filter_cells(adata,min_genes=200)
		mito_genes = adata.var_names.str.upper().str.startswith('MT-')
		adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
		adata.obs['n_counts'] = adata.X.sum(axis=1).A1
		adata = adata[adata.obs['n_genes'] < 2500, :]
		adata = adata[adata.obs['percent_mito'] < 0.05, :]
		sc.pp.normalize_per_cell(adata,key_n_counts='n_counts_all',counts_per_cell_after=ncounts)
		filter_result = sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=n_top_genes, inplace=False)
		gene_list = adata.var_names[filter_result.highly_variable]
		if hg:
			gene_list = gene_list.union(pd.read_csv("./data/marker_genes.csv").x)
		else:
		  print("No additional marker genes for non-human samples")
	else:
		gene_list=pd.read_csv(path+"genelist.csv",index_col=0).index
		
	adata = adata[:,gene_list]
	sc.pp.normalize_per_cell(adata,counts_per_cell_after=ncounts)
	# used for loss evaluation
	if filter_disp_genes:
		adata.obs['size_factors'] = (adata.obs.n_counts_all / ncounts)*(adata.obs.n_counts/ncounts)
	else:
		adata.obs['size_factors'] = adata.obs.n_counts/ncounts
	sc.pp.log1p(adata)
	if scale:
		sc.pp.scale(adata)
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
	
def normalization(adata,path,filter_disp_genes=True,hg=True,ncounts=1e6,scale=False):
	adataraw=adata.copy()
	adata=recipe_zheng(adata,path=path,filter_disp_genes=filter_disp_genes,hg=hg,ncounts=ncounts,scale=scale)
	cell_idx=adata.obs_names
	gene_idx=adata.var_names
	adata_subset = adataraw[:,gene_idx]
	adata_subset = adata_subset[cell_idx,:]
	adata.raw=adata_subset.copy()
	return adata
	
def recipe_li(adata,n_top_genes=1000):
	import desc
	sc.pp.filter_cells(adata, min_genes=200)
	sc.pp.filter_genes(adata, min_cells=3)
	mito_genes = adata.var_names.str.startswith('MT-')
	adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
	adata.obs['n_counts'] = adata.X.sum(axis=1).A1
	# sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],jitter=0.4, multi_panel=True,show=False)
	adata = adata[adata.obs['n_genes'] < 2500, :]
	adata = adata[adata.obs['percent_mito'] < 0.05, :]
	desc.normalize_per_cell(adata, counts_per_cell_after=1e4)
	desc.log1p(adata)
	adata.raw=adata
	sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, subset=True)
	adata = adata[:, adata.var['highly_variable']]
	desc.scale(adata, zero_center=True, max_value=3)
	return adata

def reference_alignment(redadata,tn=None,t=None):
	obs=redadata.obs
	obs['cid']=range(obs.shape[0])
	# col_cell_type='louvain'
	# obs['cell_type']=adata.obs[col_cell_type]
	if tn is not None:
		obstest=obs.loc[obs.replicate_id.cat.codes>=tn, ['cid','cell_type']]
	elif t is not None:
		obstest=obs.loc[obs.test==t,['cid','cell_type']]
	else:
		raise ValueError('Error: test dataset are required to be specific!')
	obstest['map_to_ref']=obstest['cell_type']
	obstest['assigned']=False
	conn_csr = redadata.uns['neighbors']['connectivities']
	for rn in obstest.index:
		# tid: test cell id; rn: barcode of cell
		tid = obs.loc[rn,'cid']
		[start, stop] = conn_csr.indptr[[tid,tid+1]]
		neighbors = conn_csr.indices[range(start,stop)]
		if len(neighbors)==0:
			continue
		dfn = obs.iloc[neighbors]
		if tn is not None:
			dfn = dfn.loc[dfn.replicate_id.cat.codes<tn]
		elif t is not None:
			dfn=dfn.loc[obs.test!=t]
		else:
			raise ValueError('Error: test dataset are required to be specific!')
			
		if dfn.empty:
			continue
		obstest.loc[rn,'map_to_ref'] = dfn.cell_type.value_counts().index[0]
		obstest.loc[rn,'assigned']=True
	return obstest

def batch_correction(obstest,redadata,tn=None,t=None,method=None):
	obstest=obstest[obstest.assigned]
	obs=redadata.obs
	if tn is not None:
		obsref=obs.loc[obs.replicate_id.cat.codes<tn, ['cid','cell_type']]
	elif t is not None:
		obsref=obs.loc[obs.test!=t,['cid','cell_type']]
	else:
		raise ValueError('Error: test dataset are required to be specific!')
		
	for ctype in obstest.map_to_ref.cat.categories.values:
		ind_test=obstest[obstest.map_to_ref==ctype].index
		ind_ref=obsref[obsref.cell_type==ctype].index
		if methon is None:
			centroidX=redadata[ind_test,:].X.mean(0).copy()
			centroidR=redadata[ind_ref,:].X.mean(0).copy()
		else:
			centroidX=np.median(redadata[ind_test,:].X)
			centroidR=np.median(redadata[ind_ref,:].X)
		centroidS=centroidR-centroidX
		redadata[ind_test,:].X=redadata[ind_test,:].X+centroidS
	return redadata
	
def nn_embedding(adata,path="./",n_neighbors=15):
	sc.pp.neighbors(adata,n_neighbors=n_neighbors)
	sc.tl.louvain(adata)
	# sc.tl.paga(adata)
	# sc.pl.paga(adata, plot=False)
	# sc.tl.umap(adata,init_pos="paga",min_dist=0.1)
	sc.tl.tsne(adata, n_jobs=20)
	adata.write(path+"ica_clusters.h5ad",compression='gzip')
	
def plotEmbedding(adata,path,color_col="louvain",ncol=10):
	# sc.pl.umap(adata,color=color_col,show=False)
	# pl.title("")
	# pl.legend(loc=3,fontsize=6,mode="expand",bbox_to_anchor=(0.0, 1.01, 1, 0.2),ncol=ncol)
	# pl.savefig(path+color_col+"_umap.png")
	# pl.close()
	sc.pl.tsne(adata,color=color_col,show=False)
	pl.title("")
	pl.legend(loc=3,fontsize=6,mode="expand",bbox_to_anchor=(0.0, 1.01, 1, 0.2),ncol=ncol)
	pl.savefig(path+color_col+"_tsne.png")
	pl.close()
def plotBatchEffect(adata,path,batch=False):
	nn_embedding(adata,path=path)
	plotEmbedding(adata,path=path,ncol=10)
	plotEmbedding(adata,path=path,color_col="cell_type_code",ncol=10)
	if batch:
		plotEmbedding(adata,path=path,color_col="batch", ncol=10)
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
		
def plotHeatmapDEG(adata,path,ntop=20,groupby="batch"):
	sc.tl.rank_genes_groups(adata, groupby=groupby, method='wilcoxon')
	result = adata.uns['rank_genes_groups']
	groups = result['names'].dtype.names
	df = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'pvals_adj']}).head(ntop)
	df.to_csv(path+"markers_"+groupby+"_"+method+".csv")
	top_ranked_genes=pd.DataFrame(adata.uns['rank_genes_groups']['names'][range(ntop)])
	gene_index = pd.Index(top_ranked_genes.values.flatten()).drop_duplicates(keep='first')
	sc.pl.heatmap(adata,gene_index, use_raw=False, swap_axes=True, show=False, show_gene_labels=True, save=False, groupby=groupby, dendrogram=False)
	pl.savefig(path+"heatmap"+groupby+"_"+method+".png")
	pl.close()

def plotQQdeg(adata,path,groupby="batch",method="wilcoxon"):
  sc.tl.rank_genes_groups(adata, groupby=groupby, method=method,n_genes=adata.shape[1],rankby_abs=True,use_raw=False)
  result = adata.uns['rank_genes_groups']
  groups = result['names'].dtype.names
  df = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'pvals_adj']})
  df.to_csv(path+"markers_"+groupby+"_"+method+".csv")
  pvals=adata.uns['rank_genes_groups']['pvals_adj']['0']+1e-260

def plotQQdeg2(adata,path,groupby="batch",method="wilcoxon"):
  for c in adata.obs.cell_type.cat.categories.values:
    adatatemp= adata[adata.obs.cell_type==c,:]
    patht=path+"pvals"+c.replace("/","_")
    plotQQdeg(adatatemp,patht)

def plotHist(adata,path):
  for c in adata.obs.cell_type.cat.categories.values:
    adatatemp= adata[adata.obs.cell_type==c,:]
    c=c.replace("/","_")
    sc.tl.rank_genes_groups(adata, groupby="batch", method="wilcoxon",n_genes=adatatemp.shape[1],rankby_abs=True,use_raw=False)
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    df = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'pvals_adj']})
    g1 = df[df['0_p']>=0.99]['0_n']
    g2 = df[df['0_p']<1e-10]['0_n']
    gl = pd.concat([g2,g1])
    import random
    sp = random.sample(list(gl.values),k=2)
    for s in sp:
      adatab1=adatatemp[adatatemp.obs.batch=='0',:]
      x1=adatab1[:,s].X
      adatab2=adatatemp[adatatemp.obs.batch=='1',:]
      x2=adatab2[:,s].X
      pl.subplot(1,2,1)
      pl.hist(x1)
      pl.subplot(1,2,2)
      pl.hist(x2)
      pl.savefig(path+"hist_"+c+"_"+s+".png")
      pl.close()
    
    
