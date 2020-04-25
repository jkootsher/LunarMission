import math
import numpy

from lib.tools.conversions import deg2rad
from lib.tools.generators import euler_232


class Kepler(object):
    ''' Kepler Elements and Calculation '''

    def __init__(self, e=0, a=1e03, raan=0, i=0, nu=0, argp=0):
        self.eccentricity = e
        self.semi_major_axis = a # km
        self.raan = deg2rad(raan)
        self.inclination = deg2rad(i)
        self.true_anomaly = deg2rad(nu)
        self.arg_at_perigee = deg2rad(argp)
    
    @classmethod
    def mean_anomaly(cls, **kwargs):
        ''' Governing dynamics for the mean anomaly '''
        e = kwargs['e']
        E = kwargs['x']
        M = kwargs['M']
        f = (E - e*math.sin(E)) - M
        return f

    @classmethod
    def mean_anomaly_df(cls, **kwargs):
        ''' Required for Newton Raphson method '''
        e = kwargs['e']
        E = kwargs['x']
        df = 1 - e*math.cos(E)
        return df

    def perifocal_to_orbital(self):
        ''' Perifocal PQW frame to orbital frame '''
        phi = (-self.true_anomaly, math.pi/2, 3*math.pi/2)
        R = euler_232(phi)
        # Each partition is 3x3 so matrix is 6x6
        # PQR2ORB = | R      ZEROS |
        #           | ZEROS  R     |
        PQR2ORB = numpy.zeros((6,6))
        PQR2ORB[:3,:3] = R
        PQR2ORB[3:,3:] = R
        return PQR2ORB

    def orbital_to_perifocal(self):
        ''' Orbital frame to perifocal PQW frame '''
        PQR2ORB = self.perifocal_to_orbital()
        return PQR2ORB.T
        