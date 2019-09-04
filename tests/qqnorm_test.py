import qqnorm
import numpy as np
a=np.matrix([[1,2,3],[4,5,6]],dtype=float)
q=qqnorm.qqnorm(a,'col')
print(q)

