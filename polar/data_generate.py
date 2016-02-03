from __future__ import absolute_import, division, print_function

import itertools
import copy
import numpy as np
from horton import *
from .fields_generate import Fields
from .function_fit import poly_fit

def finitefield_ham(ham, lf, obasis, f_order, field, xyz=None):
    """
    """

    if not xyz:
        xyz = np.zeros(3)
    
    mm = DenseTwoIndex(obasis.nbasis)
    for perm in itertools.combinations_with_replacement((0, 1, 2), f_order):
        mask = np.zeros(3, dtype=np.int64)
        for i in perm:
            mask[i] += 1
        fact = np.math.factorial(f_order)
        fact /= np.math.factorial(mask[0])*np.math.factorial(mask[1])*np.math.factorial(mask[2])
        mm_tmp = obasis.compute_multipolemoment(mask, xyz, lf)
        mm_tmp._array = mm_tmp._array*field.item(perm)
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


def model_finitefield_ham(ham, lf, obasis, olp, orb, occ_model, f_order, method='hf'):
    """
    """

    ham_backup = ham
    energys = []
    obj = Fields()
    fields = obj.fields(f_order)
    for i, field in enumerate(fields):
        ham = copy.deepcopy(ham_backup)
        ffham = finitefield_ham(ham, lf, obasis, f_order, field)
        energys.append(finitefield_energy(ffham, lf, olp, orb, occ_model, method=method))

    return poly_fit(fields, energys, 10, f_order)

