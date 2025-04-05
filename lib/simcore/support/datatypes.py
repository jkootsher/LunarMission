import math

from lib.simcore.support.variables import GLOBALS
from lib.simcore.support.exceptions import Verify
from lib.simcore.support.exceptions import VectorIndexError
from lib.simcore.support.exceptions import ReferenceFrameError

class Parameter(object):
    ''' Basic parameter type '''
    def __init__(self, value=0):
        self.ic = value
        self.ts = [value]
        return


class Vector(object):
    ''' Vector data type '''
    def __init__(self, values=[0,0,0], frame=None):
        self.n = 0
        self.frame = frame

        self.i = Parameter(values[0])
        self.j = Parameter(values[1])
        self.k = Parameter(values[2])
        return

    def __mul__(self, other=None):
        ''' Vector inner (dot) product '''
        if not isinstance(other, Vector):
            raise TypeError
        if self.frame is not other.frame:
            raise ReferenceFrameError
        
        other.n = self.n
        ii = self.i.ts[self.n]*other.i.ts[other.n]
        jj = self.j.ts[self.n]*other.j.ts[other.n]
        kk = self.k.ts[self.n]*other.k.ts[other.n]
        return ii + jj + kk

    def __pow__(self, other=None):
        ''' Vector wedge (cross) product '''
        if not isinstance(other, Vector):
            raise TypeError
        if self.frame is not other.frame:
            raise ReferenceFrameError

        cross = Vector()
        other.n = self.n
        cross.n = self.n
        cross.i = self.j.ts[self.n]*other.k.ts[other.n] - self.k.ts[self.n]*other.j.ts[other.n]
        cross.j = self.k.ts[self.n]*other.i.ts[other.n] - self.i.ts[self.n]*other.k.ts[other.n]
        cross.k = self.i.ts[self.n]*other.j.ts[other.n] - self.j.ts[self.n]*other.i.ts[other.n]
        return cross

    def __add__(self, other=None):
        ''' Vector addition '''
        if self.frame is not other.frame:
            raise ReferenceFrameError

        vector = Vector()
        vector.n = self.n

        if not isinstance(other, Vector):
            if Verify.is_numeric(other):
                vector.i.ts[vector.n] = self.i.ts[self.n] + other
                vector.j.ts[vector.n] = self.j.ts[self.n] + other
                vector.k.ts[vector.n] = self.k.ts[self.n] + other
            else:
                raise TypeError
        else:
            other.n = self.n
            vector.i.ts[vector.n] = self.i.ts[self.n] + other.i.ts[other.n]
            vector.j.ts[vector.n] = self.j.ts[self.n] + other.j.ts[other.n]
            vector.k.ts[vector.n] = self.k.ts[self.n] + other.k.ts[other.n]
        return vector

    def __sub__(self, other=None):
        ''' Vector subtraction '''
        if self.frame is not other.frame:
            raise ReferenceFrameError

        vector = Vector()
        vector.n = self.n

        if not isinstance(other, Vector):
            if Verify.is_numeric(other):
                vector.i.ts[vector.n] = self.i.ts[self.n] - other
                vector.j.ts[vector.n] = self.j.ts[self.n] - other
                vector.k.ts[vector.n] = self.k.ts[self.n] - other
            else:
                raise TypeError
        else:
            other.n = self.n
            vector.i.ts[vector.n] = self.i.ts[self.n] - other.i.ts[other.n]
            vector.j.ts[vector.n] = self.j.ts[self.n] - other.j.ts[other.n]
            vector.k.ts[vector.n] = self.k.ts[self.n] - other.k.ts[other.n]
        return vector

    def norm(self):
        ''' Calculate the vector norm '''
        n = self.n
        return math.sqrt(self.i.ts[n]**2 + self.j.ts[n]**2 + self.k.ts[n]**2)

    def at(self, dn=0):
        ''' Query the vector values at time dn '''
        return (self.i.ts[dn],self.j.ts[dn],self.k.ts[dn])

    def update(self, coordinates=[0,0,0]):
        ''' Update the vector time series '''
        if len(coordinates) > GLOBALS.CONSTANTS['VECTOR_SIZE']:
            raise VectorIndexError
        self.i.ts.append(coordinates[0])
        self.j.ts.append(coordinates[1])
        self.k.ts.append(coordinates[2])
        return


class Satellite(object):
    ''' Generic satellite class '''
    def __init__(self):
        pass