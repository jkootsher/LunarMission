import numpy

from lib.tools.conversions import mrp2dcm
from lib.tools.conversions import deg2rad
from lib.tools.math_toolbox.ode import solver
from lib.gnc.control.pointing import Pointing
from lib.gnc.control.controller import Controller


class Vehicle(object):
    ''' Default Model Parameters '''

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
    ''' Generic Orbiter Model '''

    def __init__(self, Orbit=None):
        super(Satellite, self).__init__()
        self._Controller = Controller()
        self._Pointing = Pointing(Orbit)

        # Vehicle tracking (optional)
        self._Target = None
        self._antenna_span = 0

        # Enable control
        self._control_enabled = True

        # Initial configuration
        self._Controller._moi = self.moments

    def _normalize_attitude(self, vehicle_state=None):
        ''' Check the MRPs for shadow set '''
        attitude = vehicle_state[0:3]

        # Shadow set check
        norm_s = numpy.linalg.norm(attitude)
        if norm_s > 1:
            attitude = -attitude/norm_s**2
            vehicle_state[0:3] = attitude
        return numpy.reshape(vehicle_state, (6,1))

    def _configure_antenna(self, **kwargs):
        ''' Configure the vehicle for communication '''
        # Update the antenna span
        antenna_span = self._antenna_span
        delta_angle = self._Pointing.relay.update_span(**kwargs)

        # Check if vehicle is in range
        if (delta_angle < antenna_span):
            self._Controller.set_pointing(self._Pointing.relay)
        return

    def _issue_control_command(self, n=0, **kwargs):
        ''' Issue a control command '''
        # Calculate transition points
        vector_norm = numpy.linalg.norm(kwargs['orbit'])
        e = self._Pointing.Orbit.Kepler.eccentricity
        a = self._Pointing.Orbit.Kepler.semi_major_axis
        transition = a*(1-e**2)/(1+e*numpy.cos(deg2rad(165)))

        # Primary operational modes
        transition_parameter = vector_norm < transition
        if kwargs['orbit'][0] <= 0 and transition_parameter:
            self._Controller.set_pointing(self._Pointing.sun)
        else:
            self._Controller.set_pointing(self._Pointing.nadir)

            # Verify vehicle tracking (relay mode)
            if self._Target is not None and transition_parameter:
                kwargs['nspan'] = (n,n+1)
                kwargs['local'] = self._Pointing.get_state_range(**kwargs)
                kwargs['target'] = self._Pointing.relay.track(self._Target, **kwargs)
                self._configure_antenna(**kwargs)
        
        if not (n % self._Controller._tau): # Update every second
            kwargs['control_torque'] = self._Controller.get_control_torque(**kwargs)
        return kwargs
    
    def toggle_control(self):
        ''' Toggle the controller '''
        self._control_enabled = not self._control_enabled
        return

    def configure_tracking(self, Target=None, los=42):
        ''' Configure target tracking with a specified line of sight '''
        self._Target = Target
        self._antenna_span = deg2rad(los)
        return

    def update_dynamics(self, dt=(0,1), **kwargs):
        ''' Get the state trajectory with body rates (inertial lunar frame) '''
        step_size = kwargs['delta']
        vehicle_body_state = numpy.append(self.mrp, self.rates, axis=0)

        # Initialize required conditions
        kwargs['moments'] = self.moments
        kwargs['control_torque'] = numpy.zeros((6,1))
        
        # Configure timing
        tspan = dt[-1]-dt[0]
        delta = (0,step_size)
        samples = int(1+tspan/step_size)

        # Set the controller signal decay time
        self._Controller._tau = int(1/step_size)
        
        # Vehicle state history and generation
        vehicle_state_history = vehicle_body_state

        for n in range(samples):
            kwargs['state'] = vehicle_body_state
            kwargs['orbit'] = self._Pointing.Orbit.inertial_solution(n*step_size)

            if self._control_enabled:
                kwargs = self._issue_control_command(n, **kwargs)

            # Integrate body dynamics
            [_,states] = solver.ode45v(tspan=(0,step_size), **kwargs)
            vehicle_body_state = self._normalize_attitude(states[:,-1])

            # Update vehicle state history
            kwargs['state'] = vehicle_body_state
            vehicle_state_history = numpy.append(vehicle_state_history, vehicle_body_state, axis=1)
        return vehicle_state_history