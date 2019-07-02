import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as pl
from anndata import AnnData
import sys
sys.path.insert(0,'./lib/pylib/')
import preprocess_10x as prep
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

class VAE:
	def __init__(self,input_size,path="./"):
		self.input_size=input_size
		self.path=path
		self.vae=None
		self.inputs=None
		self.outputs=None
		callbacks = []
		checkpointer = ModelCheckpoint(filepath=path+"vae_mlp_weights.h5", verbose=1, save_best_only=True)
		reduce_lr = ReduceLROnPlateau(monitor='val_loss', factor=0.1, patience=3, min_lr=0.001)
		early_stop = EarlyStopping(monitor='val_loss', patience=20)
		tensor_board = TensorBoard(log_dir=path+'logs/')
		callbacks.append(checkpointer)
		callbacks.append(reduce_lr)
		callbacks.append(early_stop)
		callbacks.append(tensor_board)
		self.callbacks = callbacks
	def sampling(self, args):
		z_mean, z_log_var = args
		batch = K.shape(z_mean)[0]
		dim = K.int_shape(z_mean)[1]
		epsilon = K.random_normal(shape=(batch, dim))
		return z_mean + K.exp(0.5 * z_log_var) * epsilon
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
		z = Lambda(self.sampling, output_shape=(32,), name='hidden_var_z')([z_mean, z_log_var])
		encoder = Model(inputs, [z_mean, z_log_var, z], name='encoder_mlp')
		# encoder.summary()
		# plot_model(encoder, to_file=path+'vae_mlp_encoder.png', show_shapes=False)
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
		# decoder.summary()
		# plot_model(decoder, to_file=path+'vae_mlp_decoder.png', show_shapes=False)
		# build vae 
		outputs = decoder(encoder(inputs)[2])
		vae = Model(inputs, outputs, name='vae_mlp')
		reconstruction_loss = mse(inputs, outputs)
		reconstruction_loss *= self.input_size
		kl_loss = 1 + z_log_var - K.square(z_mean) - K.exp(z_log_var)
		kl_loss = -0.5*K.sum(kl_loss, axis=-1)
		vae_loss = K.mean(reconstruction_loss + kl_loss)
		vae.add_loss(vae_loss)
		self.vae = vae
		self.encoder = encoder
		self.decoder = decoder
	def compile(self):
		self.vae.compile(optimizer='Adagrad')
		self.vae.summary()
		plot_model(self.vae,to_file=self.path+'vae_mlp.png')
	def train(path, batch_size=218, epochs=3000):
		self.vae.fit(adata.X, epochs=epochs, batch_size=batch_size, callbacks=net.callbacks, validation_split=0.1,shuffle=True)
		self.vae.save_weights(path+'vae_mlp_weights.h5')
		self.encoder.save_weights(path+"encoder_mlp_weights.h5")
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
	path="./result/sample/scvae_shekhar/"
	filename="./data/mouse/GSE81904/GSE81904_bipolarumicounts.h5ad"
	load_weights=True
	batch_size=512
	adata, n_batches, input_size = load_data(filename,batch_size,ref=True)
	net = VAE(input_size,path)
	net.build()
	net.compile()
	if load_weights:
		net.vae.load_weights(path+"weights_best.h5")
		net.encoder.load_weights(path+"weights_best.h5",by_name=True)
	else:
		net.train(path)
	
	[z_mean, z_log_var, z] = net.encoder.predict(adata.X)
	redadata1 = AnnData(X=z_mean, obs=adata.obs)
	path1=path+"reference_"
	prep.nn_embedding(redadata1,path=path1)
	prep.plotEmbedding(redadata1,path=path1,ncol=10)
	prep.plotEmbedding(redadata1,path=path1,color_col="cell_type_code",ncol=10)

	adata, n_batches, input_size = load_data(filename,batch_size,ref=False)
	[z_mean, z_log_var, z] = net.encoder.predict(adata.X)
	redadata2 = AnnData(X=z_mean, obs=adata.obs)
	path2=path+"test_"
	prep.nn_embedding(redadata2,path=path2)
	prep.plotEmbedding(redadata2,path=path2,ncol=10)
	prep.plotEmbedding(redadata2,path=path2,color_col="cell_type_code",ncol=10)
	
	redadata3=redadata1.concatenate(redadata2)
	path3=path+"combined_"
	prep.nn_embedding(redadata3,path=path3)
	prep.plotEmbedding(redadata3,path=path3,ncol=10)
	prep.plotEmbedding(redadata3,path=path3,color_col="cell_type_code",ncol=10)
	prep.plotEmbedding(redadata3,path=path3,color_col="batch", ncol=10)

	
