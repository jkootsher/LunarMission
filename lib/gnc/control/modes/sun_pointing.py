import numpy

from lib.gnc.control.frames import Frames
from lib.tools.conversions import mrp2dcm
from lib.tools.conversions import unit_vector


class Sun(Frames):
    ''' Sun Pointing Control Mode '''

    def __init__(self):
        super(Sun, self).__init__()
        self.fname = 'sun'

    def relative_rates(self, **kwargs):
        ''' Relative rate calculation for the sun frame '''
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
        ''' Inertial origin frame to orbiter sun pointing frame '''
        attitude = kwargs['state'][0:3,:]
        BYI2BDY = mrp2dcm(attitude)
        
        BDY2BYI = BYI2BDY.T

        # Solar array is symmetric about the body y-axis
        # and light will be incident on the body +z-axis (panel faces)
        BDY2SUN = numpy.zeros([3,3])
        BDY2SUN[0,:] = +BDY2BYI[2,:]
        BDY2SUN[1,:] = +BDY2BYI[1,:]
        BDY2SUN[2,:] = -BDY2BYI[0,:]

        # BDY2SUN*BYI2BDY*x = BYI2SUN*x
        BYI2SUN = numpy.matmul(BDY2SUN, BYI2BDY)
        return BYI2SUN

    def control_to_inertial(self, **kwargs):
        ''' Orbiter sun pointing frame to inertial origin frame '''
        BYI2SUN = self.inertial_to_control(**kwargs)
        return BYI2SUN.T