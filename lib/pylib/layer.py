from keras.engine.topology import Layer

class SliceLayer(Layer):

	def __init__(self, index, **kwargs):
		self.index = index
		super().__init__(**kwargs)

	def build(self, input_shape):
		if not isinstance(input_shape,list):
			raise ValueError('Input should be a list')
			
		super().build(input_shape)  # Be sure to call this at the end

	def call(self, x):
		assert isinstance(x, list), 'SliceLayer input is not a list'
		return x[self.index]

	def compute_output_shape(self, input_shape):
		return input_shape[self.index]
		
	def get_config(self):
		config = {'index': self.index}
		base_config = super().get_config()
		return dict(list(base_config.items()) + list(config.items()))
