import numpy as np
from keras.models import Model,Sequential
from keras.layers import Dense, Input
from keras.optimizers import Adam

class AutoEncoder():
    def __init__(self,
                 input_size,
                 ):

        self.input_size = input_size
        self.output_size = input_size

    def build(self,x_train):
        x_train_noisy = x_train + 0.01 * np.random.normal(loc=0.0, scale=1.0, size=(1,self.input_size))
        model = Sequential()
        model.add(Dense(128, activation='relu', input_dim=self.input_size))
        model.add(Dense(32, activation='relu'))
        model.add(Dense(2, activation='relu'))
        model.add(Dense(32, activation='relu'))
        model.add(Dense(128, activation='relu'))
        model.add(Dense(self.output_size, activation='relu'))

        adam = Adam(lr = 1e-4)
        model.compile(optimizer='adam', loss='mse')

        model.fit(x_train_noisy, x_train,
          	epochs=1000,
          	batch_size=128)
        score = model.evaluate(x_train_noisy, x_train, batch_size=128)
        print(score)

