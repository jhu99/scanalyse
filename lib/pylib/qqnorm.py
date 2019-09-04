from scipy.stats import norm
from scipy.stats import rankdata
import numpy as np
def cacu_theory_quantiles(a,size):
    if size > 10:
        for i in range(size):
             a[i] = (a[i]-0.5) / size
    else:
        for i in range(size):
            a[i] = (a[i]-0.375) / (size+0.25)
    for i in range(size):
        a[i] = norm.ppf(a[i])

def rank(X, method):
    for i in range(X.shape[0]):
        X[i] = rankdata(X[i], method=method)

def qqnorm(X, method='max', replace=True, axis=0):
    if replace==False:
        Y=np.zeros(X.shape,dtype=float)
        for i in range(X.shape[0]):
            for j in range(X.shape[1]):
                Y[i,j] = X[i,j]
        if axis == 0:
            rank(Y, method)
        if axis == 1:
            Y = Y.T
            rank(Y, method)
        size = Y.shape[1]
        for i in range(Y.shape[0]):
            cacu_theory_quantiles(Y[i],size)
        if axis == 1:
            Y = Y.T
        return Y
    else:
        if axis == 0:
            rank(X, method)
        if axis == 1:
            X = X.T
            rank(X, method)
        size = X.shape[1]
        for i in range(X.shape[0]):
            cacu_theory_quantiles(X[i], size)
        if axis == 1:
            X = X.T
        return X