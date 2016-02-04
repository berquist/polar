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
        list_field = []
        for j,term in enumerate(combinations_with_replacement((0,1,2),f_order)):
            list_field.append(field[j])
        row=[]
        for k in range(1,polar_order):
            for term in combinations_with_replacement(list_field,k):
                row.append(np.prod(term))
        X.append(row)
    X=np.array(X)
    return lstsq(X,bili)[0][6:30]

