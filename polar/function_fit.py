from itertools import *

import numpy as np
from numpy.linalg import lstsq
import sympy as sp
from scipy.optimize import leastsq
from math import *


def rat_fit(fields, energys, f_order, e_order):
    """
    Model some data as a rational function. It return back a list of 
    coefficient with sequnce, which is coefficient of denominator goes first
    and followed by coefficient of numerator.

    **Arguments**


    fields
        a set of fields for rational function fitting

    energys
        a set of energys for rational function fitting

    f_order
        order of field

    e_order
        expansion order. Which is a list with two elements,
        the first one is order of denominator and the second
        one if order of numerator
    """
    bili=[]
    for i in energys[1:]:
        tmp=i-energys[0]
        bili.append(tmp)

    X=[] 
    for i,field in enumerate(fields[1:]):
        row=[]
        for j in range(e_order[0]):
            for term in combinations_with_replacement(field,j+1):
                tmpp = np.prod(term)*-energys[i+1]
                row.append(tmpp)
        for k in range(e_order[1]):
            for term in combinations_with_replacement(field,k+1):
                tmpp = np.prod(term)
                row.append(tmpp)
        X.append(row)
    X=np.array(X)
    return lstsq(X,bili)[0]


def create_poly(f_order, e_order, coeff, coeff_0):
    """
    This function create a scipy polynomial functions with given coefficient

    **Arguments**


    coeff
        a list of coefficient

    coeff_0
        the first coefficent which is usually 1 for denominator E(0) for numerator
    """

    # Create the coordinate system
    coord = []
    for i in f_order:
        for term in combinations_with_replacement('xyz',i):
            coord.append(''.join(term))
    field = []
    for term in coord:
        field.append(sp.symbols('F_{}'.format(term)))
    # Create the polynomial
    poly = coeff_0
    c = 0
    for i in range(e_order):
        for term in combinations_with_replacement(field,i+1):
            tmp =1
            for bili in term:
                
                tmp *= bili
            poly += tmp*coeff[c]
            c = c+1
    return poly


def create_ratfunc(f_order, top_order, bottom_order, top_coeff, bottom_coeff, top_coeff0, bottom_coeff0):
    """
    This function create a scipy rational functions with given coefficient

    """

    # Process the coordinate system
    
    # Create the rational function
    ratfunc = create_poly(f_order, top_order, top_coeff, top_coeff0)
    ratfunc /= create_poly(f_order, bottom_order, bottom_coeff, bottom_coeff0)
    return ratfunc




def solve(fields, energys, f_order, e_order, p_order):
    """
    This function create a scipy rational functions and derivative it according 
    to the polarizability order.

    """
    # Create the coordinate system
    coord = []
    for i in f_order:
        for term in combinations_with_replacement('xyz',i):
            coord.append(''.join(term))
    field = []
    for term in coord:
        field.append(sp.symbols('F_{}'.format(term)))
    coeff = rat_fit(fields, energys, f_order, e_order)
    n0 = 0
    n1 = 0
    for i in range(e_order[0]):
        i += 1
        n0 +=factorial(i+len(field)-1)/(factorial(len(field)-1)*factorial(i))
    for j in range(e_order[1]):
        j += 1
        n1 +=factorial(j+len(field)-1)/(factorial(len(field)-1)*factorial(j))
    bottom_order = e_order[0]
    top_order = e_order[1]

    bottom_coeff = coeff[0:n0]
    top_coeff = coeff[n0:]
    bottom_coeff0 = 1
    top_coeff0 = energys[0]

    fun = create_ratfunc(f_order, top_order, bottom_order, top_coeff, bottom_coeff, top_coeff0, bottom_coeff0)
    fun_backup = fun
    polar = []
    for term in combinations_with_replacement(field, p_order):
        fun = fun_backup
        for i in term:
            fun = fun.diff(sp.symbols(str(i)))
        fun = -fun.subs({sp.symbols('F_{}'.format(i)): 0. for i in coord})
        polar.append(fun)

    return polar




