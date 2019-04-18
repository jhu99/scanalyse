import sys
sys.path.insert(0,'./lib/pylib/')
import preprocess_10x as prep
import zinb_AE as zae
import scanpy as sc
import scanpy.api
import matplotlib.pyplot as pl
import seaborn as sns
import random
sns.set(style='white', rc={'figure.figsize':(8,6), 'figure.dpi':150})
############## Read projection and plot clusters########
filepath=sys.argv[1]
adata = prep.read_10x_data(filepath+"projection.csv","10x_csv")
zae.plotCluster(adata,filepath=filepath,dm_reduction=True)
