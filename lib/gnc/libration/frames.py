import math
import numpy

from lib.simcore.support.variables import GLOBALS
# All frames (ECI, LCI, Barycentric) are relative to the
# barycentric system of Earth-Moon-Libration2. This means that
# ECI is Earth centered and inertial within this system, LCI is
# moon centered and inertial within this system, and the Barycentric
# frame is the ROTATING barycentric frame. Consider renaming frames to
# accurately reflect this. Otherwise, the language won't make sense.

class Frames(object):
    ''' Halo orbit frames class '''

    def _fixed_rotation(self, relative_frame_angle=0):
        ''' Generate the rotation matrix for a general inertial barycentric
            body centered frame to the rotating (fixed) barycentric frame '''
        FBY2BYIpos = numpy.asarray([[+math.cos(relative_frame_angle), -math.sin(relative_frame_angle), 0], \
                                    [+math.sin(relative_frame_angle), +math.cos(relative_frame_angle), 0], \
                                    [0, 0, 1]])

        FBY2BYIdot = numpy.asarray([[-math.sin(relative_frame_angle), -math.cos(relative_frame_angle), 0], \
                                    [+math.cos(relative_frame_angle), -math.sin(relative_frame_angle), 0], \
                                    [0, 0, 0]])

        ZEROS = numpy.zeros((3,3))

        # Each partition is 3x3 so matrix is 6x6
        # FBY2BYI = | FBY2BYIpos  ZEROS      |
        #           | FBY2BYIdot  FBY2BYIpos |
        FBY2BYI = numpy.zeros((6,6))
        FBY2BYI[:,:3] = numpy.concatenate((FBY2BYIpos, FBY2BYIdot), axis=0)
        FBY2BYI[:,3:] = numpy.concatenate((ZEROS, FBY2BYIpos), axis=0)
        return FBY2BYI

    def barycenter_fixed_to_barycenter_lunar_inertial(self, state_vector=None):
        ''' Rotating barycenter frame (normalized) to inertial barycenter lunar origin frame (SI units) '''
        x_barycenter = state_vector[0]
        y_barycenter = state_vector[1]

        relative_frame_angle = math.atan2(y_barycenter, x_barycenter)

        FBY2BYI = self._fixed_rotation(relative_frame_angle)

        state_vector = numpy.matmul(FBY2BYI, state_vector)
        state_vector = state_vector*GLOBALS.EARTH['TO_MOON']
        state_vector[0] = state_vector[0] - GLOBALS.EARTH['TO_MOON']
        state_vector[3:6] = state_vector[3:6]/GLOBALS.LUNAR['T_ORBIT']
        state_vector = numpy.reshape(state_vector, (6,1))
        return state_vector

    def barycenter_lunar_inertial_to_barycenter_fixed(self, state_vector=None):
        ''' Inertial barycenter lunar origin frame (SI units) to rotating barycenter frame (normalized) '''
        x_lunar = state_vector[0]
        y_lunar = state_vector[1]

        relative_frame_angle = math.atan2(y_lunar, x_lunar)

        FBY2BYI = self._fixed_rotation(relative_frame_angle)

        state_vector = numpy.matmul(FBY2BYI.T, state_vector)
        state_vector = state_vector/GLOBALS.EARTH['TO_MOON']
        state_vector[0] = state_vector[0] + 1
        state_vector[3:6] = state_vector[3:6]*GLOBALS.LUNAR['T_ORBIT']
        state_vector = numpy.reshape(state_vector, (6,1))
        return state_vector

    def barycenter_fixed_to_barycenter_earth_inertial(self, state_vector=None):
        ''' Rotating barycenter frame (normalized) to inertial barycenter Earth origin frame (SI units) '''
        x_barycenter = state_vector[0]
        y_barycenter = state_vector[1]

        relative_frame_angle = math.atan2(y_barycenter, x_barycenter)

        FBY2BYI = self._fixed_rotation(relative_frame_angle)

        state_vector = numpy.matmul(FBY2BYI, state_vector)
        state_vector = state_vector*GLOBALS.EARTH['TO_MOON']
        state_vector[3:6] = state_vector[3:6]/GLOBALS.LUNAR['T_ORBIT']
        state_vector = numpy.reshape(state_vector, (6,1))
        return state_vector

    def barycenter_earth_inertial_to_barycenter_fixed(self, state_vector=None):
        ''' Inertial barycenter Earth origin frame (SI units) to rotating barycenter frame (normalized) '''
        x_lunar = state_vector[0]
        y_lunar = state_vector[1]

        relative_frame_angle = math.atan2(y_lunar, x_lunar)

        FBY2BYI = self._fixed_rotation(relative_frame_angle)

        state_vector = numpy.matmul(FBY2BYI.T, state_vector)
        state_vector = state_vector/GLOBALS.EARTH['TO_MOON']
        state_vector[3:6] = state_vector[3:6]*GLOBALS.LUNAR['T_ORBIT']
        state_vector = numpy.reshape(state_vector, (6,1))
        return state_vector
