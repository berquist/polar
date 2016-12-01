from __future__ import absolute_import, division, print_function

from itertools import (combinations, combinations_with_replacement)
from random import uniform
from math import pi, sin, cos

import numpy as np

class Fields(object):
    def __init__(self, mol, samples=50, F=0.00005, x=2**0.5, f=None):
        '''
        This class generate fields needed for function fitting

        **Arguments:**


        samples
            number of sample points in each sphere shell

        F
            the initial field value

        x
            the interval of field
        '''
        self.samples = samples
        self.F = F
        self.x = x
        if not f:
            self.f = []
        for i in range(5):
            self.f.append(self.F * self.x**i)

    def get_l(self, mol):
        '''
        This function mesure the x, y, z distance compoent between
        all possible atom pairs and return back the biggest distance

        **Arguments:**


        mol
            object moleculer come from horton
        '''
        l = []
        for i, term in enumerate(mol.numbers):
            tmp = mol.coordinates[:, i]
            column = []
            for term in combinations(tmp, 2):
                bili = 2.6 + abs(term[0] - term[1])
                column.append(bili)
            l.append(max(column))
        return l

    def polar_coordinate(self, order):
        '''
        This function create a random point on a n dimensional hypersphere
        in polar coordinate

        **Arguments:**


        order
            field order
        '''
        phi = []
        for i in order:
            for term in combinations_with_replacement((0, 1, 2), i):
                phi.append(uniform(0, pi))
        phi = phi[0:-1]
        phi[-1] = phi[-1] * 2
        return phi

    def sphere_list(self, order, r):
        '''
        This funciton return back a n dimensional hyper sphere field list

        **Arguments:**


        order
            field order

        r
            a list of radius of several sphere layers
            [F0, F0*x, F0*x**2, F0*x**3 ...]
        '''
        field = []
        phi = self.polar_coordinate(order)
        for i, angle in enumerate(phi):
            if i == 0:
                x = r
                x = x * cos(phi[i])
                field.append(x)
            elif i > 0 and i < len(phi) - 1:
                x = r
                for j in range(i):
                    x = x * sin(phi[j])
                x = x * cos(phi[i])
                field.append(x)
            elif i == len(phi) - 1:
                x = r
                for j in range(i):
                    x = x * sin(phi[j])
                x = x * cos(phi[i])
                field.append(x)
                x = r
                for j in range(i + 1):
                    x = x * sin(phi[j])
                field.append(x)

        return field

    def ellipse_list(self, order, r, mol):
        '''
        This funciton return back a n dimensional hyper ellipse field list

        **Arguments:**


        order
            field order

        r
            a list of radius of several sphere layers
            [F0, F0*x, F0*x**2, F0*x**3 ...]

        mol
            object moleculer come from horton
        '''
        l = self.get_l(mol)
        field = self.sphere_list(order, r)
        c = 0
        for i in order:
            for term in combinations_with_replacement((0, 1, 2), i):
                tmp = 1
                for j in term:
                    tmp *= l[j]
                field[c] *= tmp
                c += 1
        return field

    def fields(self, order, mol):
        '''
        This function return back a set of fields with certain field
        pattern

        **Arguments:**


        order
            field order

        mol
            object moleculer come from horton
        '''
        n = 0
        for i in order:
            n += (i+1)*(i+2)/2
        fields = [np.zeros(n)]
        for r in self.f:
            for i in range(self.samples):
                tmp = self.ellipse_list(order, r, mol)
                fields.append(np.array(tmp))
        return fields
