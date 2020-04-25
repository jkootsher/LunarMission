import math
import numpy

from lib.tools.generators import euler_313
from lib.simcore.support.variables import GLOBALS


class Frames(object):
    ''' Lunar orbit frames class '''

    def __init__(self, Kepler=None):    
        p0 = Kepler.arg_at_perigee
        p1 = Kepler.inclination
        p2 = Kepler.raan
        
        # Rotation angle set
        self.phi = (p0, p1, p2)

    def barycenter_lunar_inertial_to_perifocal(self):
        ''' Inertial barycenter lunar origin frame to perifocal PQW frame '''
        R = euler_313(self.phi).T
        # Each partition is 3x3 so matrix is 6x6
        # BYI2PQW = | R      ZEROS |
        #           | ZEROS  R     |
        BYI2PQW = numpy.zeros((6,6))
        BYI2PQW[:3,:3] = R
        BYI2PQW[3:,3:] = R
        return BYI2PQW

    def perifocal_to_barycenter_lunar_inertial(self):
        ''' Perifocal PQW frame to inertial barycenter lunar origin frame '''
        R = euler_313(self.phi)
        # Each partition is 3x3 so matrix is 6x6
        # PQW2BYI = | R      ZEROS |
        #           | ZEROS  R     |
        PQW2BYI = numpy.zeros((6,6))
        PQW2BYI[:3,:3] = R
        PQW2BYI[3:,3:] = R
        return PQW2BYI

    def barycenter_lunar_inertial_to_barycenter_earth_inertial(self, state_vector=None):
        ''' Inertial barycenter lunar origin frame to inertial barycenter Earth origin frame '''
        offset = GLOBALS.EARTH['TO_MOON']
        state_vector[0,0] = state_vector[0,0] + offset
        return state_vector

    def barycenter_earth_inertial_to_barycenter_lunar_inertial(self, state_vector=None):
        ''' Inertial barycenter Earth origin frame to inertial barycenter lunar origin frame '''
        offset = GLOBALS.EARTH['TO_MOON']
        state_vector[0,0] = state_vector[0,0] - offset
        return state_vector