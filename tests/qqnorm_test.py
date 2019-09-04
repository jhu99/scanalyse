import sys
sys.path.insert(0,'./lib/pylib/')
import qqnorm
import numpy as np
b=np.array([[1,2,3],[4,5,6]],ndmin=2,dtype=float)
q=qqnorm.qqnorm(b,'max',True,0)
print(b)

b=np.array([[1,2,3],[4,5,6]],ndmin=2,dtype=float)
q=qqnorm.qqnorm(b,'max',False,0)
print(b)
print(q)

