from __future__ import absolute_import, division, print_function

import itertools

import numpy as np
from horton import *


def finitefield_ham(ham, lf, obasis, order_1, order_2, field_1, field_2, xyz=None):
    """
    """

    if not xyz:
        xyz = np.zeros(3)


    mm = DenseTwoIndex(obasis.nbasis)
    for perm in itertools.combinations_with_replacement((0, 1, 2), order_1):
        mask_1 = np.zeros(3, dtype=np.int64)
        for i in perm:
            mask_1[i] += 1
        fact = np.math.factorial(order_1)
        fact /= np.math.factorial(mask_1[0])*np.math.factorial(mask_1[1])*np.math.factorial(mask_1[2])
        mm_tmp = obasis.compute_multipolemoment(mask_1, xyz, lf)
        mm_tmp._array = mm_tmp._array*field_1.item(perm)
        mm._array += mm_tmp._array*fact
    
    if order_1==order_2:
        pass
    else:
        for perm in itertools.combinations_with_replacement((0, 1, 2), order_2):
            mask_2 = np.zeros(3, dtype=np.int64)
            for i in perm:
                mask_2[i] += 1
            fact = np.math.factorial(order_2)
            fact /= np.math.factorial(mask_2[0])*np.math.factorial(mask_2[1])*np.math.factorial(mask_2[2])
            mm_tmp = obasis.compute_multipolemoment(mask_2, xyz, lf)
            mm_tmp._array = mm_tmp._array*field_2.item(perm)
            mm._array += mm_tmp._array*fact

    ham.terms.append(RTwoIndexTerm(mm,'mm'))

    return ham


def finitefield_energy(ham, lf, olp, orb, occ_model, method='hf'): 
    """
    """

    if method is 'hf':
        scf_solver = PlainSCFSolver(1e-6)
        scf_solver(ham, lf, olp, occ_model, orb)
    else:
        raise Exception

    return ham.compute_energy()
