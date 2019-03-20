import sys
sys.path.insert(0,'../lib/pylib/')
from preprocess_10x import getNormAnnData
import AE

if(sys.argv[2]=='log'):
    adata,size_factor = getNormAnnData(sys.argv[1],sys.argv[2])
    AE.train_zinb_model(adata,size_factor)
    AE.prediction_zinb_middle(adata,size_factor,sys.argv[3])
else:
    adata= getNormAnnData(sys.argv[1],sys.argv[2])
    AE.train_zip_model(adata)
    AE.prediction_zip_middle(adata,sys.argv[3])
