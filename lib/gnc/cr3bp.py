import numpy

from lib.simcore.support.variables import GLOBALS

class LCR3BP(object):
    ''' Parameter Class for the CR3BP (L2 only for now) '''

    def __init__(self, T=2.4499, h=1e-4, sbdy=None, lbdy=None, cycles=1):
        self.position = numpy.asarray([1.0874,0.0,0.2020])
        self.velocity = numpy.asarray([0.0,-0.2054,0.0])
        self.cycles = cycles
        self.period = T
        self.step = h

        # Set the orbital bodies (MU)
        self.sbdy = sbdy
        self.lbdy = lbdy

    def change_initial_conditions(self, r=None, v=None, T=0):
        ''' Change the initial orbital parameters '''
        self.position = numpy.asarray(r)
        self.velocity = numpy.asarray(v)
        self.period = T
        return

    def period_in_seconds(self):
        ''' Get the period in seconds '''
        return GLOBALS.LUNAR['T_ORBIT']*self.period

    def total_samples(self):
        ''' Get the total sample count '''
        return int(1+self.period/self.step)