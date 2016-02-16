from __future__ import absolute_import, division, print_function
from horton import *
import numpy as np
from itertools import *
from math import *
from random import *

class Fields(object):
    def __init__(self, mol, samples=50, F=0.005, x=2**0.5 ,f=None):
    
        '''
        This class generate fields needed for function fitting

        **Arguments:**

        order
            A tuple with two field order

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
        if f is None:
            self.f=[]
        for i in range(5):
            self.f.append(self.F*self.x**i)

    def get_l(self, mol):
        l = []
        for i,term in enumerate(mol.numbers):
            tmp = mol.coordinates[:,i]
            column = []
            for term in combinations(tmp , 2):
                bili = abs(term[0]-term[1])
                if bili==0.:
                    bili=0.48
                column.append(bili)
            l.append(max(column))
        l = [1,1,1]
        return l

    def polar_coordinate(self, order):
        phi=[]
        for i in order:
            for term in combinations_with_replacement((0,1,2), i):
                phi.append(uniform(0,pi))
        phi=phi[0:-1]
        phi[-1]=phi[-1]*2
        return phi


    def sphere_list(self, order, r):
        '''
        This funciton return back a filed list 
        '''
        field = []
        phi = self.polar_coordinate(order)
        for i,angle in enumerate(phi):
            if i==0:
                x=r
                x=x*cos(phi[i])
                field.append(x)
            elif i>0 and i<len(phi)-1:
                x=r
                for j in range(i):
                    x=x*sin(phi[j])
                x=x*cos(phi[i])
                field.append(x)
            elif i==len(phi)-1:
                x=r
                for j in range(i):
                    x=x*sin(phi[j])
                x=x*cos(phi[i])
                field.append(x)
                x=r
                for j in range(i+1):
                    x=x*sin(phi[j])
                field.append(x)

        return field

    def ellipse_list(self, order, r, mol):
        l = self.get_l(mol)
        field = self.sphere_list(order,r)
        c = 0
        for i in order:
            for term in combinations_with_replacement((0,1,2),i):
                tmp = 1
                for j in term:
                    tmp *= l[j]
                field[c] *=tmp
                c += 1
        return field

    def fields(self,order,mol):
        '''
        '''
        n = 0
        for i in order:
            n += (i+1)*(i+2)/2
        fields=[np.zeros(n)]
        for r in self.f:
            for i in range(self.samples):
                tmp = self.ellipse_list(order, r, mol)
                fields.append(np.array(tmp))
        return fields
