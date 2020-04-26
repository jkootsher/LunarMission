#!/usr/bin/env python

from lib.gnc.cr3bp import LCR3BP
from lib.gnc.kepler import Kepler
from lib.models.orbiter import Satellite
from lib.models.lunar_relay import Target
from lib.simcore.support.variables import GLOBALS
from lib.gnc.lunar.dynamics import Orbit as Lunar
from lib.gnc.libration.dynamics import Orbit as L2

class Simulator(object):
    ''' Simulator '''

    def run(self):
        elements = Kepler(e=0.6, a=3000, i=56.2, argp=90, raan=0)
        Lunar_Orbit_A = Lunar(elements, phase=0)
        Lunar_Orbit_B = Lunar(elements, phase=120)
        Lunar_Orbit_C = Lunar(elements, phase=240)

        Orbiter_A = Satellite(Lunar_Orbit_A)
        Orbiter_B = Satellite(Lunar_Orbit_B)
        Orbiter_C = Satellite(Lunar_Orbit_C)