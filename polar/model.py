#from __future__ import absolute_import, division, print_function

from copy import deepcopy
import numpy as np
from finitefield.horton import finitefield_ham, finitefield_energy
from ratfunc.solve import model
from random import *
from math import *

def fields_creator(f_order,samples):
    F=0.00005
    x=2**0.5
    f=[]
    for i in range(5):
        f.append(F*x**i)
    fields=[np.zeros((3,)*f_order)]
    shape=[]
    for i in range(f_order):
        shape.append(3)

#################################
    for r in f:
        for n in range(samples):
            field=np.zeros((3,)*f_order)
            for i in range(3**f_order):
                if i==0:
                    x=r*cos(uniform(0,1*pi))
                    field.append(x)
                elif i==3**f_order:
                    x=r
                    for j in range(i-1):
                        x=x*sin(uniform(0,1*pi))
                    x=x*sin(uniform(0,2*pi))
                    field.append(x)
                else:
                    x=r
                    for j in range(i-1):
                        x=x*sin(uniform(0,1*pi))
                    x=x*cos(uniform(0,2*pi))
                    field.append(x)
            fields.append(np.array(field).reshape(shape)) 
    return fields



def model_finitefield_ham(order, ham, lf, obasis, olp, orb, occ_model, method='hf'):
    """
    """

    ham_backup = ham

    samples = 0
    for i in range(1,order):
        samples = samples + (i+8)*(i+7)*(i+6)*(i+5)*(i+4)*(i+3)*(i+2)*(i+1)/40320
    energys = []
    fields=fields_creator(2,50)
    fields=fields[0:(samples+1)]

    
    for i, field in enumerate(fields):
        ham = deepcopy(ham_backup)
        ffham = finitefield_ham(ham, lf, obasis, field=field)
        energys.append(finitefield_energy(ffham, lf, olp, orb, occ_model, method=method))
    return model(order, fields, energys)

