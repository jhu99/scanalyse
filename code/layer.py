from keras.engine.topology import Layer
from keras.layers import Lambda, Dense
import tensorflow as tf

class SliceLayer(Layer):
    def __init__(self, index, **kwargs):
        self.index = index
        super().__init__(**kwargs)

    def build(self, input_shape):
        if not isinstance(input_shape, list):
            raise ValueError('Input should be a list')

        super().build(input_shape)

    def call(self, x):
        assert isinstance(x, list), 'SliceLayer input is not a list'
        return x[self.index]

    def compute_output_shape(self, input_shape):
        return input_shape[self.index]

ColwiseMultLayer = Lambda(lambda l: l[0]*tf.reshape(l[1], (-1,1)))