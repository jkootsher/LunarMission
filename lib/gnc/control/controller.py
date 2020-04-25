import numpy

from lib.gnc.control import modes
from lib.tools.conversions import mrp2dcm
from lib.tools.conversions import dcm2mrp
from lib.gnc.control.modes.stub_pointing import Stub


class Controller(object):
    ''' Simple PD controller '''

    def __init__(self, tau=10):
        self.tau = tau              # Decay time
        self.moi = None             # Principal moments

        # Controls
        self.set_pointing()

    def set_pointing(self, mode=Stub()):
        self._mode = mode
        self._frame = mode.inertial_to_control
        return

    def get_control_torque(self, **kwargs):
        ''' Obtain the necessary control torque '''
        state_estimate = self._estimate_current_state(**kwargs)
        
        kp = numpy.linalg.norm(self.moi/self.tau, numpy.inf)
        ku = 2*numpy.linalg.norm(self.moi/self.tau**2, numpy.inf)

        MOI_INV = numpy.linalg.inv(self.moi)
        state_estimate = ku*state_estimate[0:3] + kp*state_estimate[3:6]
        state_estimate = numpy.reshape(state_estimate, (3,1))
        return -numpy.matmul(MOI_INV, state_estimate)

    def _estimate_current_state(self, **kwargs):
        ''' Estimate the relative state for the given control frame '''
        vehicle_states = kwargs['state']

        # Ascertain state information
        mrp_attitude = numpy.reshape(vehicle_states[0:3], (3,1))
        rates_vector = numpy.reshape(vehicle_states[3:6], (3,1))

        # Construct the inertial to control frame
        BYI2BDY = mrp2dcm(mrp_attitude)
        BYI2REF = self._frame(**kwargs)

        # Body to control frame
        BDY2REF = numpy.matmul(BYI2REF, BYI2BDY.T)
        
        # Body update
        mrp_attitude = dcm2mrp(BDY2REF.T)
        rates_vector = numpy.matmul(BYI2BDY.T, rates_vector) - \
                       self._mode.relative_rates(**kwargs)
        rates_vector = numpy.matmul(BYI2BDY, rates_vector)

        vehicle_bdy_update = numpy.empty((6,1))
        vehicle_bdy_update[0:3,:] = mrp_attitude
        vehicle_bdy_update[3:6,:] = rates_vector
        return vehicle_bdy_update
    