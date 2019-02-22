import numpy as np
from sklearn.preprocessing import StandardScaler
import tensorflow as tf
from tensorflow.examples.tutorials.mnist import input_data


class AddtiveGussianNoiseAutoencoder(object):
    def __init__(self, n_input, n_hidden1, n_hidden2, scale, optimizer, transfer_function=tf.nn.softplus):
        self.n_input = n_input
        self.n_hidden1 = n_hidden1
        self.n_hidden2 = n_hidden2
        self.transfer = transfer_function
        self.optimizer = optimizer
        self.scale = tf.placeholder(tf.float32)
        self.training_scale = scale
        network_weight = self._initialize_weights()
        self.weights = network_weight

        self.x = tf.placeholder(tf.float32, [None, self.n_input])
        self.hidden1 =self.transfer(tf.add(tf.matmul(self.x + scale * tf.random_normal([n_input]), self.weights['w1']),
                                           self.weights['b1']))
        self.hidden2=self.transfer(tf.add(tf.matmul(self.hidden1 , self.weights['w2']),
                                           self.weights['b2']))
        self.reconstruction1 = self.transfer(tf.add(tf.matmul(self.hidden2, self.weights['w3']),
                                                   self.weights['b3']))
        self.reconstruction2 = self.transfer(tf.add(tf.matmul(self.reconstruction1, self.weights['w4']),
                                                   self.weights['b4']))
        self.cost = 0.5 * tf.reduce_sum(tf.pow(tf.subtract(self.reconstruction2, self.x), 2.0))
        self.optimizer = optimizer.minimize(self.cost)

        init = tf.global_variables_initializer()
        self.sess = tf.Session()
        self.sess.run(init)

    def _initialize_weights(self):
        all_weights = dict()
        all_weights['w1'] = tf.Variable(xavier_init(self.n_input, self.n_hidden1))
        all_weights['b1'] = tf.Variable(tf.zeros([self.n_hidden1], dtype=tf.float32))
        all_weights['w2'] = tf.Variable(xavier_init(self.n_hidden1, self.n_hidden2))
        all_weights['b2'] = tf.Variable(tf.zeros([self.n_hidden2], dtype=tf.float32))
        all_weights['w3'] = tf.Variable(xavier_init(self.n_hidden2, self.n_hidden1))
        all_weights['b3'] = tf.Variable(tf.zeros([self.n_hidden1], dtype=tf.float32))
        all_weights['w4'] = tf.Variable(xavier_init(self.n_hidden1, self.n_input))
        all_weights['b4'] = tf.Variable(tf.zeros([self.n_input], dtype=tf.float32))
        return all_weights

    def partial_fit(self, X):
        cost, opt = self.sess.run((self.cost, self.optimizer),
                                  feed_dict={self.x: X,
                                             self.scale: self.training_scale})
        return cost

    def calc_total_cost(self, X):
        return self.sess.run((self.cost, self.optimizer),
                             feed_dict={self.x: X,
                                        self.scale: self.training_scale})

    def generate(self, hidden=None):
        if hidden is None:
            hidden = np.random.normal(size=self.weights['w1'].shape[0])
        return self.sess.run(self.reconstruction, feed_dict={self.hidden: hidden})

    def reconstruct(self, X):
        return self.sess.run(self.reconstruction, feed_dict={self.x: X, self.scale: self.training_scale})

    def get_weights(self):
        return self.sess.run(self.weights['w1'])

    def get_biases(self):
        return self.sess.run(self.weights['w2'])


def standard_scale(X_train, X_test):
    scaler = StandardScaler()
    scaler.fit(X_train)
    X_train = scaler.transform(X_train)
    scaler.fit(X_test)
    X_train = scaler.transform(X_test)
    return X_train, X_test


def get_random_block_from_data(data, batch_size):
    start_index = np.random.randint(0, n_sample-batch_size)
    return data[start_index: (start_index+batch_size)]

def xavier_init(fan_in, fan_out, constant=1):
    low = -constant * np.sqrt(6.0 / (fan_in + fan_out))
    high = constant * np.sqrt(6.0 / (fan_in + fan_out))
    return tf.random_uniform((fan_in, fan_out), minval=low, maxval=high,
                             dtype=tf.float32)


if __name__ == '__main__':
    mnist = input_data.read_data_sets('MNIST_data', one_hot=True)
    x_train, x_test = standard_scale(mnist.train.images, mnist.train.images)
    n_sample = int(mnist.train.num_examples)

    batch_size = 128
    display_step = 1
    train_epoch = 200

    autoencoder = AddtiveGussianNoiseAutoencoder(
        n_input=784,
        n_hidden1 = 64,
        n_hidden2=32,
        transfer_function=tf.nn.softplus,
        optimizer=tf.train.AdamOptimizer(learning_rate=0.001),
        scale  = 0.1
    )

    for epoch in range(train_epoch):
        avg_cost = 0
        total_batch = int(n_sample/batch_size)
        for i in range(total_batch):
            batch_xs = get_random_block_from_data(x_train, batch_size)
            cost = autoencoder.partial_fit(batch_xs)
            avg_cost += cost / n_sample * batch_size

        if epoch % display_step == 0:
            print("Epoch: %04d" % (epoch + 1), "cost={:.9f}".format(avg_cost))
    print("Total cost:" + str(autoencoder.calc_total_cost(x_test)))