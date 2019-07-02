import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as pl
from anndata import AnnData
import sys
sys.path.insert(0,'./lib/pylib/')
import preprocess_10x as prep
from vae import VAE
from keras.layers import Input, Dense, Activation, BatchNormalization, Dropout, Lambda
from keras.models import Model, Sequential, load_model
from keras.callbacks import TensorBoard, ModelCheckpoint, EarlyStopping, ReduceLROnPlateau
from keras import backend as K
from keras import optimizers, regularizers
from keras.utils import plot_model
from keras.losses import mse, binary_crossentropy
import seaborn as sns

matplotlib.use('Agg')
sns.set(style='white', rc={'figure.figsize':(8,6), 'figure.dpi':150})

def sampling(args):
	z_mean, z_log_var = args
	batch = K.shape(z_mean)[0]
	dim = K.int_shape(z_mean)[1]
	epsilon = K.random_normal(shape=(batch, dim))
	return z_mean + K.exp(0.5 * z_log_var) * epsilon
def vae_loss_function(y_true, y_pred):
	reconstruction_loss = mse(y_true, y_pred)
	vae_loss_ = K.mean(reconstruction_loss)
	return vae_loss_		
class AAE:
	def __init__(self,input_size,batch_size,path="./"):
		self.input_size=input_size
		self.path=path
		self.batch_size=batch_size
		self.aae=None
		self.inputs=None
		self.outputs=None
	def build_discriminator(self):
		model = Sequential()
		model.add(Dense(16, input_dim=32, activation="relu", name="disc_layer1"))
		model.add(Dense(8, activation="relu", name="disc_layer2"))
		model.add(Dense(1, activation="sigmoid", name="disc_layer3"))
		return model
	def build(self):
		# build encoder
		Relu="relu"
		inputs = Input(shape=(self.input_size,), name='encoder_input')
		x = Dense(128)(inputs)
		x = BatchNormalization(center=True,scale=False)(x)
		x = Activation(Relu)(x)
		x = Dense(64)(x)
		x = BatchNormalization(center=True,scale=False)(x)
		x = Activation(Relu)(x)
		z_mean = Dense(32, name="encoder_mean")(x)
		z_log_var = Dense(32, name="encoder_log_var")(x)
		z = Lambda(sampling, output_shape=(32,), name='hidden_var_z')([z_mean, z_log_var])
		encoder = Model(inputs, [z_mean, z_log_var, z], name='encoder_mlp')
		# build decoder
		latent_inputs = Input(shape=(32,), name='z_sampling')
		x = Dense(64)(latent_inputs)
		x = BatchNormalization(center=True,scale=False)(x)
		x = Activation(Relu)(x)
		x = Dense(128)(x)
		x = BatchNormalization(center=True,scale=False)(x)
		x = Activation(Relu)(x)
		outputs = Dense(self.input_size, activation='softplus')(x)
		decoder = Model(latent_inputs,outputs, name='decoder_mlp')
		# build varational autoencoder
		outputs = decoder(encoder(inputs)[2])
		vae = Model(inputs, outputs, name='vae_mlp')
		# build discriminator and aae
		disc = self.build_discriminator()
		disc.trainable=False
		outputs = disc(encoder(inputs)[2])
		aae = Model(inputs,outputs,name='aae_mlp')
		# build generator of a prior distribution 
		gen = VAE(self.input_size,self.path)
		gen.build()
		gen.encoder.load_weights(self.path+"encoder_mlp_weights.h5")
		# build aae
		self.vae = vae
		self.aae = aae
		self.encoder = encoder
		self.decoder = decoder
		self.disc = disc
		self.gen = gen
	def generate_samples(self,radata,batch_size=256):
		nindex = np.random.choice(radata.obs.index,size=batch_size)
		samples = radata[nindex,:]
		[r_mean,r_log_var,r_batch] = self.gen.encoder.predict(samples.X)
		return [r_mean,r_log_var,r_batch]
	def compile(self):
		self.vae.compile(optimizer='Adagrad', loss=vae_loss_function)
		self.aae.compile(optimizer='sgd', loss='binary_crossentropy')
		self.disc.trainable=True
		self.disc.compile(optimizer='sgd', loss='binary_crossentropy')
		self.gen.compile()
		self.vae.summary()
		self.aae.summary()
		self.disc.summary()
	def predict(self, adata, radata):
		[r_mean,r_log_var,r_batch] = self.gen.encoder.predict(radata.X)
		[z_mean, z_log_var, z] = self.encoder.predict(adata.X)
		r_adata = AnnData(X=r_mean, obs=radata.obs)
		z_adata = AnnData(X=z_mean, obs=adata.obs)
		c_adata = r_adata.concatenate(z_adata)
		prep.plotBatchEffect(r_adata, path=self.path+"reference_")
		prep.plotBatchEffect(z_adata, path=self.path+"test_")
		prep.plotBatchEffect(c_adata, path=self.path+"combined_", batch=True)
		
	def load_weights(self):
		self.encoder.load_weights(self.path+"encoder_aae_weights.h5")
	def train_aae(self, adata, radata, n_batches, path, batch_size=256, n_epochs=3000):
		for epoch in range(n_epochs):
			# nindex=np.random.permutation(adata.obs.index)
			# adata = adata[nindex,:]
			for index in range(n_batches):
				# for each mini-batch, update autoencoder by minimazing the reconstruction error.
				x_batch = adata.X[index*batch_size:(index+1)*batch_size]
				x_batch = x_batch.toarray()
				vae_loss = self.vae.train_on_batch(x_batch,x_batch)
				# update discriminator
				self.disc.trainable = True
				[z_mean,z_log_var,z_batch] = self.encoder.predict(x_batch)
				[r_mean,r_log_var,r_batch] = self.generate_samples(radata,batch_size=batch_size)
				combined_batch = np.concatenate((r_batch, z_batch))
				y_true = [1] * batch_size + [0] * batch_size
				d_loss = self.disc.train_on_batch(combined_batch,y_true)
				# update encoder
				self.disc.trainable = False
				aae_loss = self.aae.train_on_batch(x_batch, [1] * batch_size)
				if epoch%199 == 0:
					print("epoch %d batch %d, vae_loss, d_loss, aae_loss : %f, %f, %f" % (epoch, index, vae_loss, d_loss, aae_loss))
			if epoch%199 == 0:
				self.encoder.save_weights(self.path+"encoder_aae_weights.h5")
				self.decoder.save_weights(self.path+"decoder_aae_weights.h5")
def load_data(filename,batch_size,ref=True):
	adata = prep.read_10x_data(filename,format_type="10x_h5ad")
	adata.obs['cell_type_code']=adata.obs['cell_type'].cat.codes
	adata.obs['cell_type_code']=adata.obs['cell_type_code'].astype("category")
	if ref:
		adata = adata[adata.obs.replicate_id.cat.codes<4,:]
	else:
		adata = adata[adata.obs.replicate_id.cat.codes>=4,:]
	adata = prep.normalization(adata,path=path,hg=False,filter_disp_genes=ref)
	n_batches = int(adata.shape[0]/batch_size)
	input_size = int(adata.shape[1])
	nindex=np.random.permutation(adata.obs.index)
	adata = adata[nindex,:]
	return adata, n_batches, input_size

if __name__ == '__main__':
	path="./result/sample/scaae_shekhar/"
	filename="./data/mouse/GSE81904/GSE81904_bipolarumicounts.h5ad"
	load_weights=False
	batch_size=256
	n_epochs=30000
	radata, r_n_batches, r_input_size = load_data(filename,batch_size,ref=True)
	adata, n_batches, input_size = load_data(filename,batch_size,ref=False)
	net = AAE(input_size, batch_size = batch_size, path=path)
	net.build()
	net.compile()
	if load_weights:
		net.load_weights()
	else:
		net.train_aae(adata=adata, radata=radata, n_epochs=n_epochs, 	n_batches=n_batches, batch_size=batch_size, path=path)
	net.predict(adata=adata, radata=radata)
