from itertools import *

import numpy as np
from numpy.linalg import lstsq
import sympy as sp
from scipy.optimize import leastsq
############################################################
def poly_fit(fields, energys, polar_order, f_order):
    """
    Model some data as a rational function.

    """
    bili=[]
    for i in energys[1:]:
        tmp=i-energys[0]
        bili.append(tmp)

    X=[] 
    for i,field in enumerate(fields[1:]):
        row=[]
        for k in range(1,polar_order):
            for term in combinations_with_replacement(field,k):
                row.append(np.prod(term))
        X.append(row)
    X=np.array(X)
    return lstsq(X,bili)[0]

def rat_fit(fields, energys, polar_order, f_order):
    """
    Model some data as a rational function.

    """
    bili=[]
    for i in energys[1:]:
        tmp=i-energys[0]
        bili.append(tmp)

    X=[] 
    for i,field in enumerate(fields[1:]):
        row=[]
        for k in range(1,polar_order[0]):
            for term in combinations_with_replacement(field,k):
                tmpp = np.prod(term)*-bili[i]
                row.append(tmpp)
        for k in range(1,polar_order[1]):
            for term in combinations_with_replacement(field,k):
                tmpp = np.prod(term)
                row.append(tmpp)
        X.append(row)
    X=np.array(X)
    return lstsq(X,bili)[0]

