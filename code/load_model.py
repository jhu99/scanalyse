import numpy as np
import os
import pandas as pd
import scanpy.api as sc
import tensorflow as tf
import keras.optimizers as opt
from keras.callbacks import TensorBoard, ModelCheckpoint, EarlyStopping, ReduceLROnPlateau
from keras import backend as K
from sklearn.model_selection import train_test_split
from model import ZINBAutoencoder
from getAnnData import getAnnData


def load_weight(data_path, weight_path, result_path):
    #adata = getAnnData(data_path)
    adata = getAnnData(data_path)
    # delete gene and cell with all 0 value
    sc.pp.filter_genes(adata, min_counts=1)
    sc.pp.filter_cells(adata, min_counts=1)
    adata.raw = adata.copy()

    # calculate size factors
    # normalization
    sc.pp.normalize_per_cell(adata)
    adata.obs['size_factors'] = adata.obs.n_counts / np.median(adata.obs.n_counts)
    # log transfer and normalization
    sc.pp.log1p(adata)
    sc.pp.scale(adata)

    output_size = adata.n_vars
    input_size = adata.n_vars
    hidden_size = [64, 32, 64]
    hidden_dropout = 0

    net = ZINBAutoencoder(input_size=input_size,
                          output_size=output_size,
                          hidden_size=hidden_size,
                          l2_coef=0.0,
                          l1_coef=0.0,
                          l2_enc_coef=0.0,
                          l1_enc_coef=0.0,
                          ridge=0.0,
                          hidden_dropout=0.0,
                          input_dropout=0.0,
                          batchnorm=True,
                          activation='relu',
                          init='glorot_uniform',
                          debug=False,
                          file_path=result_path)
    net.build()
    net.load_weights(weight_path)
    net.predict(adata, mode='full', return_info=True)
    net.write(adata, mode='full')


