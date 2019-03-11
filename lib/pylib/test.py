import sys
from preprocess_10x import getNormAnnData
import AE

if sys.argv[2] == "log":
	adata,size_factor = getNormAnnData(sys.argv[1],sys.argv[2])
	AE.train_zinb_model(adata,size_factor)
	AE.prediction_zinb(adata,size_factor)

else:
	adata = getNormAnnData(sys.argv[1],sys.argv[2])
	AE.train_zip_model(adata)
	AE.prediction_zip(adata)
