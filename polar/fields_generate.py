from __future__ import absolute_import, division, print_function
from horton import *
import numpy as np
from itertools import *
from math import *
from random import *

class Fields(object):
    def __init__(self, samples=50, F=0.00005, x=2**0.25, f=None):
    
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

    def polar_coordinate(self, order):
        phi=[]
        for term in combinations_with_replacement((0,1,2), order[0]):
            phi.append(uniform(0,pi))
        if order[0]==order[1]:
            pass
        else:
            for term in combinations_with_replacement((0,1,2), order[1]):
                phi.append(uniform(0,pi))
        phi=phi[0:-1]
        phi[-1]=phi[-1]*2
        return phi


    def field_list(self, order, r):
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

    def fields_list(self,order):
        '''
        '''
        if order[0]==order[1]:
            n = (order[0]+1)*(order[0]+2)/2
        else:
            n = (order[0]+1)*(order[0]+2)/2+(order[1]+1)*(order[1]+2)/2
        fields=[np.zeros(n)]
        for r in self.f:
            for i in range(self.samples):
                tmp = self.field_list(order,r)
                fields.append(np.array(tmp))
        return fields

