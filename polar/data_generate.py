from __future__ import absolute_import, division, print_function

import itertools
import copy
import numpy as np
from horton import *
from .fields_generate import Fields
from .function_fit import poly_fit,rat_fit, solve

def finitefield_ham(ham, lf, obasis, f_order, field, xyz=None):
    """
    This function return back hamiltonian with multipole moment term in external field

    **Arguments**
    

    lf
        LinalgFactory object from horton

    obasis
        Gaussian basis set object from horton

    f_order
        order of field

    field
        field

    xyz
        center of multipole moment, it's [0. 0. 0.] by default
    """

    if not xyz:
        xyz = np.zeros(3)
    
    mm = DenseTwoIndex(obasis.nbasis)
    for j, perm in enumerate(itertools.combinations_with_replacement((0, 1, 2), f_order)):
        mask = np.zeros(3, dtype=np.int64)
        for i in perm:
            mask[i] += 1
        mm_tmp = obasis.compute_multipolemoment(mask, xyz, lf)
        mm_tmp._array = mm_tmp._array*field[j]
        mm._array += mm_tmp._array
    ham.terms.append(RTwoIndexTerm(mm,'mm'))

    return ham


def finitefield_energy(ham, lf, olp, orb, occ_model, method='hf'): 
    """
    This function return back energy with multipole moment term in external field

    **Arguments**

    
    olp
        overlap matrix from horton

    orb
        alpha orbitals from horton, be carefully it also named exp_alpha. Make sure the 
        varible name keep consistant when you use it.

    occ_model
        horton object to decide how to occupy the orbitals
    """

    if method is 'hf':
        scf_solver = PlainSCFSolver(1e-6)
        scf_solver(ham, lf, olp, occ_model, orb)
    else:
        raise Exception

    return ham.compute_energy()


def model_finitefield_ham(ham, lf, obasis, olp, orb, occ_model, mol, f_order, p_order, method='hf'):
    """
    This function return back a set of fields and the energys for function fitting

    **Arguments**

    f_order
        a list include field order

    p_order
        polarizability order

    """
    ham_backup = ham
    energys = []
    obj = Fields(mol)
    fields = obj.fields(f_order, mol)
    for i, field in enumerate(fields):
        n0 = 0
        n1 = 0
        for i in f_order:
            n1 +=(i+1)*(i+2)/2
            ham = copy.deepcopy(ham_backup)
            field_i=field[n0:n1]
            ffham = finitefield_ham(ham, lf, obasis, i, field_i)
        n0 = n1
        energys.append(finitefield_energy(ffham, lf, olp, orb, occ_model, method=method))
    return solve(fields, energys, f_order, [3,4], p_order)

