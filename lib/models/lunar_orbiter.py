import numpy

from lib.tools.math_toolbox.ode import solver
from lib.tools.conversions import mrp2dcm

from lib.gnc.control.pointing import Pointing
from lib.gnc.control.controller import Controller
from lib.tools.conversions import deg2rad


class Vehicle(object):
    ''' Default Lunar Model Parameters '''

    def __init__(self):
        self.mrp = self._default_mrps()
        self.rates = self._default_rates()
        self.moments = self._default_mois()

    def _default_mrps(self):
        ''' Default MRP set (body relative to inertial) '''
        mrp_set = numpy.array([0.3, -0.4, 0.5])
        return numpy.reshape(mrp_set, (3,1))

    def _default_rates(self):
        ''' Default rates (body relative to inertial) '''
        r_i = deg2rad(+1.00)
        r_j = deg2rad(+1.75)
        r_k = deg2rad(-2.20)
        rates = numpy.array([r_i, r_j, r_k])
        return numpy.reshape(rates, (3,1))

    def _default_mois(self):
        ''' Default moments of inertia (principal frame) '''
        return numpy.array([[100, 0, 0], [0, 25, 0], [0, 0, 75]])


class Satellite(Vehicle):
    ''' Lunar Model Class '''

    def __init__(self, Orbit=None):
        super(Satellite, self).__init__()
        self.Pointing = Pointing(Orbit)
        self.Controller = Controller()

        # Initial configuration
        self.Controller.moi = self.moments

    def _normalize_attitude(self, vehicle_state=None):
        ''' Check the MRPs for shadow set '''
        attitude = vehicle_state[0:3]

        # Shadow set check
        norm_s = numpy.linalg.norm(attitude)
        if norm_s > 1:
            attitude = -attitude/norm_s**2
            vehicle_state[0:3] = attitude
        return numpy.reshape(vehicle_state, (6,1))
        
    def update_dynamics(self, dt=(0,1), **kwargs):
        ''' Get the state trajectory with body rates (inertial lunar frame) '''
        step_size = kwargs['delta']
        vehicle_body_state = numpy.append(self.mrp, self.rates, axis=0)

        # Initialize required conditions
        kwargs['moments'] = self.moments
        
        tspan = dt[-1]-dt[0]
        samples = int(tspan/step_size)
        vehicle_state_history = vehicle_body_state

        # Set the controller signal decay time
        tau = int(1/step_size)
        self.Controller.tau = tau
        
        for n in range(samples):
            kwargs['state'] = vehicle_body_state
            kwargs['orbit'] = self.Pointing.Orbit.inertial_solution(n*step_size)

            if kwargs['orbit'][0] > 0:
                self.Controller.set_pointing(self.Pointing.nadir)
            else:
                self.Controller.set_pointing(self.Pointing.sun)
            
            if not (n % tau): # Update every second
                kwargs['control_torque'] = self.Controller.get_control_torque(**kwargs)

            [_,states] = solver.ode45v(tspan=(0,step_size), **kwargs)
            vehicle_body_state = self._normalize_attitude(states[:,-1])

            # Update vehicle state history
            kwargs['state'] = vehicle_body_state
            vehicle_state_history = numpy.append(vehicle_state_history, vehicle_body_state, axis=1)
        return vehicle_state_history