import numpy

from lib.gnc.control.modes.sun_pointing import Sun
from lib.gnc.control.modes.nadir_pointing import Nadir
from lib.gnc.control.modes.relay_pointing import Communication


class Pointing(object):
    ''' Pointing Dynamics Base Class '''

    def __init__(self, Orbit=None):
        self.Orbit = Orbit

        # Pointing Modes
        self.sun = Sun()
        self.nadir = Nadir()
        self.relay = Communication()

    def get_state_range(self, **kwargs):
        ''' Returns the range of states dn '''
        dn = kwargs['nspan']
        samples = 1+(dn[-1]-dn[0])
        state_subset = numpy.empty((6,0))

        for n in range(samples):
            state = self.Orbit.inertial_solution(n*kwargs['delta'])
            state_subset = numpy.append(state_subset, state, axis=1)
        return state_subset