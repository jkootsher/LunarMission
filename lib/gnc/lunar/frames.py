import math
import numpy

from lib.simcore.support.variables import GLOBALS


class Frame(object):
    ''' Polar orbit frames class '''

    def __init__(self, Kepler=None):
        self._KOE = Kepler

    def _get_rmat(self):
        ''' Generate the rotation matrix for Lunar ECI to PQW '''
        ZEROS = numpy.zeros((3,3))

        # Euler 3-1-3
        THETA = self._KOE.arg_at_perigee
        PHI = self._KOE.inclination
        PSI = self._KOE.raan

        r11 = math.cos(PSI)*math.cos(THETA) - math.sin(PSI)*math.cos(PHI)*math.sin(THETA)
        r12 = -math.cos(PSI)*math.sin(THETA) - math.sin(PSI)*math.cos(PHI)*math.cos(THETA)
        r13 = math.sin(PSI)*math.sin(PHI)
        r21 = math.sin(PSI)*math.cos(THETA) + math.cos(PSI)*math.cos(PHI)*math.sin(THETA)
        r22 = -math.sin(PSI)*math.sin(THETA) + math.cos(PSI)*math.cos(PHI)*math.cos(THETA)
        r23 = -math.cos(PSI)*math.sin(PHI)
        r31 = math.sin(PHI)*math.sin(THETA)
        r32 = math.sin(PHI)*math.cos(THETA)
        r33 = math.cos(PHI)

        R = numpy.asarray([[r11, r12, r13], [r21, r22, r23], [r31, r32, r33]])

        # Each partition is 3x3 so matrix is 6x6
        # LCI2PQW = | R      ZEROS |
        #           | ZEROS  R     |
        LCI2PQW = numpy.zeros((6,6))
        LCI2PQW[:,:3] = numpy.concatenate((R, ZEROS), axis=0)
        LCI2PQW[:,3:] = numpy.concatenate((ZEROS, R), axis=0)
        return LCI2PQW

    def pqw2lci(self, state_vector=None):
        ''' Perifocal PQW frame to lunar centered inertial frame '''
        LCI2PQW = self._get_rmat()
        state_vector = numpy.matmul(LCI2PQW.T, state_vector)
        return numpy.reshape(state_vector, (6,1))

    def lci2pqw(self, state_vector=None):
        ''' Lunar centered inertial frame to perifocal PQW frame '''
        LCI2PQW = self._get_rmat()
        state_vector = numpy.matmul(LCI2PQW, state_vector)
        return numpy.reshape(state_vector, (6,1))

    def lci2eci(self, state_vector=None):
        ''' Lunar centered inertial frame to Earth centered inertial frame '''
        offset = GLOBALS.EARTH['TO_MOON']
        state_vector[0,0] = state_vector[0,0] + offset
        return state_vector

    def eci2lci(self, state_vector=None):
        ''' Earth centered inertial frame to lunar centered inertial frame '''
        offset = GLOBALS.EARTH['TO_MOON']
        state_vector[0,0] = state_vector[0,0] - offset
        return state_vector