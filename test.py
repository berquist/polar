#!/usr/bin/env python
#JSON {"lot": "RHF/3-21G",
#JSON  "scf": "PlainSCFSolver",
#JSON  "linalg": "DenseLinalgFactory",
#JSON  "difficulty": 1,
#JSON  "description": "Basic RHF example with dense matrices"}


from horton import *
import numpy as np
import sympy as sp
from polar.data_generate import *
from polar.function_fit import *

# Hartree-Fock calculation
# ------------------------

# Load the coordinates from file.
# Use the XYZ file from HORTON's test data directory.
fn_xyz = context.get_fn('test/H2O.xyz')
mol = IOData.from_file(fn_xyz)
# Create a Gaussian basis set
obasis = get_gobasis(mol.coordinates, mol.numbers, 'sto-3g')
lf =CholeskyLinalgFactory(obasis.nbasis)
# Create a linalg factory
# Compute Gaussian integrals
olp = obasis.compute_overlap(lf)
kin = obasis.compute_kinetic(lf)
na = obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, lf)
er = obasis.compute_electron_repulsion(lf)



# Create alpha orbitals
orb = lf.create_expansion()

# Initial guess
guess_core_hamiltonian(olp, kin, na, orb)

# Construct the restricted HF effective Hamiltonian
external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
terms = [
    RTwoIndexTerm(kin, 'kin'),
    RDirectTerm(er, 'hartree'),
    RExchangeTerm(er, 'x_hf'),
    RTwoIndexTerm(na, 'ne'),
]
ham = REffHam(terms, external)

# Decide how to occupy the orbitals (5 alpha electrons)
occ_model = AufbauOccModel(5)
#print finitefield_ham(ham, lf, obasis, [1,1] , [np.array([1,1,1])])
print model_finitefield_ham(ham, lf, obasis, olp, orb, occ_model, 2)



