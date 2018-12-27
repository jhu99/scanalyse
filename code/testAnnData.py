from getAnnData import getAnnData
import train

path = '../data/ica_bone_marrow_h5.h5'
adata = getAnnData(path)
#print(adata)
train.train_model(adata)
