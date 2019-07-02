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
from keras.layers import Input, Dense, Activation, BatchNormalization, Dropout, Lambda
from keras.models import Model, Sequential, load_model
from keras.callbacks import TensorBoard, ModelCheckpoint, EarlyStopping, ReduceLROnPlateau
from keras import backend as K
from keras import optimizers, regularizers
import seaborn as sns

matplotlib.use('Agg')
sns.set(style='white', rc={'figure.figsize':(8,6), 'figure.dpi':150})

class GAN:
	def __init__(self,input_size):
		self.generator=None
		self.discriminator=None
		self.gan=None
		self.input_size=input_size
	def build_generator(self):
		model = Sequential()
		model.add(Dense(self.input_size, input_dim=self.input_size, activity_regularizer=regularizers.l1_l2(0.01,0.01), name="gen_layer1"))
		model.add(Activation("relu",name="act_layer1"))
		model.add(Dense(self.input_size, activity_regularizer=regularizers.l1_l2(0.01,0.01), name="gen_layer2"))
		model.add(Activation("relu",name="act_layer2"))
		model.add(Dense(self.input_size, activity_regularizer=regularizers.l1_l2(0.01,0.01), name="gen_layer3"))
		model.add(Activation("relu",name="act_layer3"))
		model.summary()
		return model
	def build_discriminator(self):
		model = Sequential()
		model.add(Dense(16, input_dim=self.input_size, activation="relu", name="disc_layer1"))
		model.add(Dense(8, activation="relu", name="disc_layer2"))
		model.add(Dense(1, activation="sigmoid", name="disc_layer3"))
		model.summary()
		return model
	def generator_containing_discriminator(self):
		model = Sequential()
		model.add(self.g)
		self.d.trainable =False
		model.add(self.d)
		model.summary()
		return model
	def build(self):
		self.g = self.build_generator()
		self.d = self.build_discriminator()
		self.gan = self.generator_containing_discriminator()
	def load_saved_model(self, path):
		self.g = load_model(path+"generator_best_model.h5")
		self.d = load_model(path+"discriminator_best_model.h5")
		self.gan = load_model(path+"gan_best_model.h5")
	def compile(self):
		d_optim = optimizers.SGD(lr=0.001, momentum=0.1, nesterov=True)
		g_optim = optimizers.SGD(lr=0.001, momentum=0.1, nesterov=True)
		# g_optim = optimizers.Adam(lr=0.001, beta_1=0.9, beta_2=0.999, epsilon=None, decay=0.0, amsgrad=False)
		self.g.compile(loss='binary_crossentropy', optimizer="sgd")
		self.gan.compile(loss='binary_crossentropy',optimizer=g_optim)
		self.d.trainable = True
		self.d.compile(loss='binary_crossentropy',optimizer=d_optim)
def load_data(filename,batch_size):
	import sys
	sys.path.insert(0,'./lib/pylib/')
	import preprocess_10x as prep
	adata = prep.read_10x_data(filename,format_type="10x_h5ad")
	# adata = adata[adata.obs.replicate_id.cat.codes<4,:]
	n_batches = int(adata.shape[0]/batch_size)
	input_size = int(adata.shape[1])
	nindex=np.random.permutation(adata.obs.index)
	adata = adata[nindex,:]
	return adata, n_batches, input_size
def train(path, n_epoch=100, batch_size=218):
	filename=path+"reference_ica_clusters.h5ad"
	adata, n_batches, input_size = load_data(filename,batch_size)
	net = GAN(input_size)
	net.build()
	net.compile()
	np.random.seed(2019)
	for epoch in range(n_epoch):
		print("Epoch is", epoch)
		nindex=np.random.permutation(adata.obs.index)
		adata = adata[nindex,:]
		for index in range(n_batches):
			noise = np.random.normal(0, 1, size=(batch_size,input_size))
			data_batch = adata.X[index*batch_size:(index+1)*batch_size]
			generated_batch = net.g.predict(noise)
			combined_batch = np.concatenate((data_batch,generated_batch))
			y = [1] * batch_size + [0] * batch_size
			# train on each batch
			net.d.trainable=True
			d_loss = net.d.train_on_batch(combined_batch, y)
			noise = np.random.normal(0, 1, (batch_size, input_size))
			net.d.trainable = False
			g_loss = net.gan.train_on_batch(noise, [1] * batch_size)
			if epoch%19 == 0:
				print("batch %d g_loss : %f" % (index, g_loss))
				print("batch %d d_loss : %f" % (index, d_loss))
		# save network model and weights
		if epoch%19 == 0:
			net.g.save_weights(path+'generator_best_weights.h5')
			net.d.save_weights(path+'discriminator_best_weight.h5')
def generate(path,n_samples=512, input_size=32):
	net = GAN(input_size)
	net.build()
	net.compile()
	net.g.load_weights(path+"generator_best_weights.h5")
	# generate random samples
	noise = np.random.normal(0, 1, (n_samples, input_size))
	generated_samples = net.g.predict(noise)
	adata = AnnData(X=generated_samples)
	adata.write(path+"generated_samples.h5")
	return adata

if __name__ == '__main__':
	path="./result/sample/sckms_shekhar/"
	train(path, n_epoch=30000)
