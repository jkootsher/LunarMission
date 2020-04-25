import numpy

from lib.tools.conversions import unit_vector


class Frames(object):
    ''' Local Navigation Frames '''

    def inertial_to_lvlh(self, state_vector=None):
        ''' A genertic inertial frame to the LVLH (Hill) frame '''
        r_vector = numpy.reshape(state_vector[0:3], (3,1))
        v_vector = numpy.reshape(state_vector[3:6], (3,1))
        w_vector = numpy.cross(r_vector, v_vector, axis=0)
        
        BYI2LVH = numpy.zeros((3,3))
        BYI2LVH[2,:] = unit_vector(r_vector.T)
        BYI2LVH[1,:] = unit_vector(w_vector.T)
        BYI2LVH[0,:] = numpy.cross(BYI2LVH[2,:], BYI2LVH[0,:])
        return BYI2LVH

    def lvlh_to_inertial(self, state_vector=None):
        ''' LVLH (Hill) frame to a generic inertial frame '''
        BYI2LVH = self.lvlh_to_inertial(state_vector)
        return BYI2LVH.T