from keras.layers import Input, Dense, Activation, BatchNormalization, Dropout, Lambda
from keras.models import Model, load_model
from keras.regularizers import l1_l2
from keras.callbacks import TensorBoard, ModelCheckpoint, EarlyStopping, ReduceLROnPlateau
from keras import backend as K

import tensorflow as tf
import pandas as pd
import csv
import numpy as np

import sys
sys.path.insert(0,'./lib/pylib/')
from layer import SliceLayer
from loss import ZINB

from keras.utils.generic_utils import get_custom_objects

# ------------------------------------------------------------
# needs to be defined as activation class otherwise error
# AttributeError: 'Activation' object has no attribute '__name__'	
class MeanAct(Activation):
	
	def __init__(self, activation, **kwargs):
		super(MeanAct, self).__init__(activation, **kwargs)
		self.__name__ = 'meanact'

def meanact(x):
	return (tf.clip_by_value(K.exp(x), 1e-5, 1e6))

class DispAct(Activation):
	
	def __init__(self, activation, **kwargs):
		super(DispAct, self).__init__(activation, **kwargs)
		self.__name__ = 'dispact'

def dispact(x):
	return (tf.clip_by_value(tf.nn.softplus(x), 1e-4, 1e4))


get_custom_objects().update({'dispact': DispAct(dispact)})
get_custom_objects().update({'meanact': MeanAct(meanact)})
get_custom_objects().update({'SliceLayer': SliceLayer})
get_custom_objects().update({'tf': tf})

def setSession():
	session_conf = tf.ConfigProto(intra_op_parallelism_threads=1, inter_op_parallelism_threads=1)
	sess = tf.Session(graph=tf.get_default_graph(), config=session_conf)
	K.set_session(sess)
	
class ZINB_AutoEncoder:
	def __init__(self):
		self.model = None
		self.loss = None
		self.encoder_model = None
		callbacks = []
		checkpointer = ModelCheckpoint(filepath='./result/ica_all/weights_best.h5', verbose=1, save_best_only=True)
		reduce_lr = ReduceLROnPlateau(monitor='val_loss', factor=0.1, patience=2, min_lr=0.001)
		early_stop = EarlyStopping(monitor='val_loss', patience=4)
		tensor_board = TensorBoard(log_dir='./result/ica_all/logs')
		callbacks.append(checkpointer)
		callbacks.append(reduce_lr)
		callbacks.append(early_stop)
		callbacks.append(tensor_board)
		self.callbacks = callbacks
	
	def build(self, input_size):
		inputs = Input(shape=(input_size,), name="counts")
		sf = Input(shape=(1,), name='size_factors')
		Relu='relu'
		
		# Construct network layers
		x = Dense(128,  kernel_regularizer=l1_l2(l1=0.,l2=0.), name="encoder_layer_1")(inputs)
		x = BatchNormalization(center=True,scale=False)(x)
		x = Activation(Relu, name="activation_el_1")(x)
		x = Dense(64, kernel_regularizer=l1_l2(l1=0.,l2=0.), name="encoder_layer_2")(x)
		x = BatchNormalization(center=True,scale=False)(x)
		x = Activation(Relu, name="activation_el_2")(x)
		x = Dense(32, kernel_regularizer=l1_l2(l1=0.,l2=0.), name="center_layer")(x)
		x = BatchNormalization(center=True,scale=False)(x)
		c = Activation(Relu, name="activation_cl")(x)
		x = Dense(64, kernel_regularizer=l1_l2(l1=0.,l2=0.), name="decoder_layer_1")(c)
		x = BatchNormalization(center=True,scale=False)(x)
		x = Activation(Relu, name="activation_dl_1")(x)
		x = Dense(128, kernel_regularizer=l1_l2(l1=0.,l2=0.), name="decoder_layer_2")(x)
		x = BatchNormalization(center=True,scale=False)(x)
		x = Activation(Relu, name="activation_dl_2")(x)
		pi= Dense(input_size, kernel_regularizer=l1_l2(l1=0.,l2=0.), activation='sigmoid', name="pi_layer")(x)
		dp= Dense(input_size, kernel_regularizer=l1_l2(l1=0.,l2=0.), activation='dispact', name="dispersion_layer")(x)
		mu= Dense(input_size, kernel_regularizer=l1_l2(l1=0.,l2=0.), activation='meanact', name="mean_layer")(x)
		ColwiseMultLayer = Lambda(lambda l: l[0]*tf.reshape(l[1], (-1,1)))
		outputs = ColwiseMultLayer([mu,sf])
		outputs = SliceLayer(0, name='slice')([outputs,dp])
		
		# Define loss function and callbacks strategies
		zinb = ZINB(pi, theta=dp, ridge_lambda=0, debug=False)
		self.loss=zinb.loss
		
		# Define models
		self.model = Model(inputs=[inputs,sf],outputs=outputs)
	
	def compile(self):
		self.model.compile(optimizer='adam', loss=self.loss)
		
	def predict(self, adata, mode='latent'):
		self.encoder_model = Model(inputs=self.model.input, outputs=self.model.get_layer('activation_cl').output)
		if mode == 'latent':
			adata.obsm['X_m'] = self.encoder_model.predict({'counts':adata.X,'size_factors':adata.obs.size_factors}, batch_size=128)
		else:
			raise ValueError("mode Error")
		return adata
		
	def write(self, adata, filepath):
		#filename = './result/ica_all/latent.csv'
		colnames = None
		rownames = adata.obs_names
		pd.DataFrame(adata.obsm['X_m'], index=rownames, columns=colnames).to_csv(filepath,sep=',',index=(rownames is not None), header=(colnames is not None), float_format='%.4f')
		
		
def train_zinb_model(adata):
	# Input and output data
	input_data={'counts':adata.X,'size_factors':adata.obs.size_factors}
	output_label=adata.raw.X
	
	# build and train the model
	net = ZINB_AutoEncoder()
	net.build(input_size=adata.n_vars)
	net.compile()
	net.model.summary()
	losses = net.model.fit(input_data, output_label, callbacks=net.callbacks, epochs=30, batch_size=128, shuffle=True, validation_split=0.1, verbose=2)
	net.model.save("./result/ica_all/model_best.h5")
	
def prediction(adata, filepath):
	setSession()
	net = ZINB_AutoEncoder()
	net.build(adata.n_vars)
	net.model.load_weights("./result/ica_all/weights_best.h5")
	net.model.summary()
	net.predict(adata)
	net.write(adata, filepath)

def plotCluster(adata,h5adfile,umapfile,tsnefile):
	import scanpy as sc
	import matplotlib.pyplot as pl
	sc.pp.neighbors(adata)
	sc.tl.louvain(adata)
	sc.tl.leiden(adata)
	sc.tl.umap(adata)
	sc.tl.tsne(adata,n_pcs=20)
	adata.write(h5adfile,compression='gzip')
	sc.pl.umap(adata,color=['louvain','leiden'],show=False)
	pl.savefig(umapfile)
	pl.close()
	sc.pl.tsne(adata,color=['louvain','leiden'],show=False)
	pl.savefig(tsnefile)
	pl.close()
	return adata
