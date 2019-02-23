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
from getAnnData import getAnnData, getAnnData_10x_h5, getAnnData_10x_mtx, pre_process_input_data
from train import train

def load_weight(input_file, weight_file_1, weight_file_2, optimizer, hidden_size, filtered, gene_file,output_path, mode, format_type="10x_h5"):
    adata = pre_process_input_data(gene_file,input_file,filtered,format_type)
    adata.raw = adata.copy()

    # calculate size factors
    # normalization
    sc.pp.normalize_per_cell(adata)
    adata.obs['size_factors'] = adata.obs.n_counts / np.median(adata.obs.n_counts)
    # log transfer and normalization
    sc.pp.log1p(adata)
    sc.pp.scale(adata,zero_center=False)

    output_size = adata.n_vars
    input_size = adata.n_vars
    #hidden_size = [64, 32, 64]
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
                          file_path=output_path)
    net.build()
    net.load_weights(weight_file_1)
    
    losses = train(adata, net, optimizer=optimizer, weight_file=weight_file_2, output_dir=output_path)
    net.predict(adata, mode=mode, return_info=True)
    net.write(adata, mode= mode)
    #net.write(adata, mode='full')


