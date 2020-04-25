import numpy

from lib.gnc.control.frames import Frames
from lib.tools.conversions import unit_vector


class Nadir(Frames):
    ''' Nadir Pointing Control Mode '''

    def __init__(self):
        super(Nadir, self).__init__()
        self.fname = 'nadir'

    def relative_rates(self, **kwargs):
        ''' Relative rate calculation for the nadir frame '''
        state_vector = kwargs['orbit']
        r_vector = numpy.reshape(state_vector[0:3], (3,1))
        v_vector = numpy.reshape(state_vector[3:6], (3,1))

        # Transform to LVLH (Hill) frame
        BYI2LVH = self.inertial_to_lvlh(state_vector)
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

        BYI2LVH = self.inertial_to_lvlh(state_vector)

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