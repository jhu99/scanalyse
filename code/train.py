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
from getAnnData import getAnnData, getAnnData_10x_h5, getAnnData_10x_mtx

def train(adata, network, output_dir=None, optimizer='rmsprop', learning_rate=0.001,
          epochs=300, reduce_lr=10, output_subset=None, use_raw_as_output=True,
          early_stop=15, batch_size=32, clip_grad=5., save_weights=True,
          validation_split=0.1, tensorboard=False, verbose=True, threads=None,
          **kwds):

    #K.set_session(tf.Session(config=tf.ConfigProto(intra_op_parallelism_threads=threads, inter_op_parallelism_threads=threads)))
    #K.clear_session()
    model = network.model
    loss = network.loss
    if output_dir is not None:
        os.makedirs(output_dir, exist_ok=True)

    optimizer = opt.rmsprop(lr=learning_rate, clipvalue=clip_grad)
    model.compile(loss=loss, optimizer=optimizer)

    # Callbacks
    callbacks = []

    checkpointer = ModelCheckpoint(filepath="%s/weights.hdf5" % output_dir,
                                   verbose=verbose,
                                   save_weights_only=save_weights,
                                   save_best_only=True)
    callbacks.append(checkpointer)

    if reduce_lr:
        lr_cb = ReduceLROnPlateau(monitor='val_loss', patience=reduce_lr, verbose=verbose)
        callbacks.append(lr_cb)
    if early_stop:
        es_cb = EarlyStopping(monitor='val_loss', patience=early_stop, verbose=verbose)
        callbacks.append(es_cb)
    if tensorboard:
        tb_log_dir = os.path.join(output_dir, 'tb')
        tb_cb = TensorBoard(log_dir=tb_log_dir, histogram_freq=1, write_grads=True)
        callbacks.append(tb_cb)

    if verbose: model.summary()

    inputs = {'count': adata.X, 'size_factors': adata.obs.size_factors}


    output = adata.raw.X

    loss = model.fit(inputs, output,
                     epochs=epochs,
                     batch_size=batch_size,
                     shuffle=True,
                     callbacks=callbacks,
                     validation_split=validation_split,
                     verbose=verbose,
                     **kwds)

    return loss

def train_model(input_file,output_path,format_type='10x_h5'):
    #K.set_session(tf.Session())
    # load data
    if format_type=="10x_h5":
        adata = getAnnData_10x_h5(input_file)
    elif format_type == "10x_mtx":
        adata = getAnnData_10x_mtx(input_file)
    else :
        raise ValueError('`format` needs to be \'10x_h5\' or \'10x_mtx\'')
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
    hidden_size = [64,32,64]
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
    losses = train(adata, net, output_dir=output_path)

    #net.predict(adata, mode='full', return_info=True)
    #net.write(adata, "./result", mode='full')

