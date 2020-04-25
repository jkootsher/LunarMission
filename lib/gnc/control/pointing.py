import numpy

from lib.gnc.control.modes.sun_pointing import Sun
from lib.gnc.control.modes.nadir_pointing import Nadir


class Pointing(object):
    ''' Pointing Dynamics Base Class '''

    def __init__(self, Orbit=None):
        self.Orbit = Orbit

        # Pointing Modes
        self.sun = Sun()
        self.nadir = Nadir()