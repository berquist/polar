from itertools import *

import numpy as np
from numpy.linalg import lstsq
import sympy as sp
from scipy.optimize import leastsq
############################################################
def model(order, fields, energys):
    """
    Model some data as a rational function.

    """
    bili=[]
    for i in energys[1:]:
        tmp=i-energys[0]
        bili.append(tmp)

    X=[] 
    for i in fields[1:]:
        view=i.ravel()
        row=[]
        for j in range(1,order):
            for term in combinations_with_replacement(view,j):
                row.append(np.prod(term))
        X.append(row)
    X=np.array(X)
    return lstsq(X,bili)[0][9:18]

