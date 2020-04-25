import numpy

from lib.tools.conversions import unit_vector


class Nadir(object):
    ''' Nadir Pointing Control Mode '''

    def __init__(self):
        self.fname = 'nadir'

    def _inertial_to_lvlh(self, state_vector=None):
        ''' General inertial to LVLH (Hill) '''
        r_vector = numpy.reshape(state_vector[0:3], (3,1))
        v_vector = numpy.reshape(state_vector[3:6], (3,1))
        w_vector = numpy.cross(r_vector, v_vector, axis=0)
        
        BYI2LVH = numpy.zeros((3,3))
        BYI2LVH[2,:] = unit_vector(r_vector.T)
        BYI2LVH[1,:] = unit_vector(w_vector.T)
        BYI2LVH[0,:] = numpy.cross(BYI2LVH[2,:], BYI2LVH[0,:])
        return BYI2LVH

    def relative_rates(self, **kwargs):
        ''' Relative rate calculation for the nadir frame '''
        state_vector = kwargs['orbit']
        r_vector = numpy.reshape(state_vector[0:3], (3,1))
        v_vector = numpy.reshape(state_vector[3:6], (3,1))

        # Transform to LVLH (Hill) frame
        BYI2LVH = self._inertial_to_lvlh(state_vector)
        r_vector = numpy.matmul(BYI2LVH, r_vector)
        v_vector = numpy.matmul(BYI2LVH, v_vector)

        # Unit vectors
        r_vector = unit_vector(r_vector)
        v_vector = unit_vector(v_vector)
        w_vector = numpy.cross(r_vector, v_vector, axis=0)

        rates = numpy.matmul(BYI2LVH.T, w_vector)
        return rates

    def inertial_to_control(self, **kwargs):
        ''' Inertial origin frame to orbiter nadir pointing frame '''
        state_vector = kwargs['orbit']

        BYI2LVH = self._inertial_to_lvlh(state_vector)

        # Antenna normally facing away from nadir
        # To point toward nadir, simply reflect x
        BYI2NAD = numpy.zeros([3,3])
        BYI2NAD[0,:] = -BYI2LVH[0,:]
        BYI2NAD[1,:] = +BYI2LVH[1,:]
        BYI2NAD[2,:] = numpy.cross(BYI2NAD[0,:], BYI2NAD[1,:])
        return BYI2NAD

    def control_to_inertial(self, **kwargs):
        ''' Orbiter nadir pointing frame to inertial origin frame '''
        BYI2NAD = self.inertial_to_control(**kwargs)
        return BYI2NAD.T