import math
import lib.tools.conversions as convert

class Kepler(object):
    ''' Kepler Elements and Calculation '''

    def __init__(self, e=0, a=1e03, raan=0, i=0, nu=0, argp=0):
        self.eccentricity = e
        self.semi_major_axis = a # km
        self.raan = convert.deg2rad(raan)
        self.inclination = convert.deg2rad(i)
        self.true_anomaly = convert.deg2rad(nu)
        self.arg_at_perigee = convert.deg2rad(argp)
    
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