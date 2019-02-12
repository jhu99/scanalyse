import numpy as np
from keras.models import Model,Sequential
from keras.layers import Dense, Input
from keras.optimizers import optm

class AutoEncoder():
    def __init__(self,
                 input_size,
                 hidden_size,
                 output_size
                 ):

        self.input_size = input_size
        self.hidden_size = hidden_size
        self.output_size = input_size

    def build(self,x_train):
        #x_train_noisy = x_train + 0.01 * np.random.normal(loc=0.0, scale=1.0, size=(1,self.input_size))
        model = Sequential()
        model.add(Dense(units=128, activation='relu', input_dim=self.input_size))
        model.add(Dense(units=32, activation='relu'))
        model.add(Dense(units=128, activation='relu'))
        model.add(Dense(self.output_size, activation='softmax'))

        adam = optm.Adam(lr = le-4, clipvalue=5)
        model.compile(optimizer=adam, loss='sparse_categorical_crossentropy')

        model.fit(x_train, x_train,
            epochs=10,
            batch_size=128, validation_split=0.1, verbose=verbose)
        
        #score = model.evaluate(x=x_train, x_train, batch_size=128)
        print(score)

