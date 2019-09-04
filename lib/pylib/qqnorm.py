from scipy.stats import norm
from scipy.stats import rankdata

def cacu_theory_quantiles(a):
    size=a.shape[1]
    if size > 10:
        for i in range(size):
             a[0,i] = (a[0,i]-0.5) / size
    else:
        for i in range(size):
            a[0, i] = (a[0, i]-0.375) / (size+0.25)
    for i in range(size):
        a[0,i] = norm.ppf(a[0,i])
    return a

def rank(X):
    for i in range(X.shape[0]):
        X[i]=rankdata(X[i], method='max')
    return X

def qqnorm(X, cacuby):
    if cacuby=='row':
        X=rank(X)
        for i in range(X.shape[0]):
            X[i,:]=cacu_theory_quantiles(X[i,:])
    else:
        X=X.T
        X=rank(X)
        for i in range(X.shape[0]):
            X[i,:] = cacu_theory_quantiles(X[i])
    return  X