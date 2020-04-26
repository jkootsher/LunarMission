import numpy


class Target(object):
    ''' Generic Relay Target Orbiter '''

    def __init__(self, Orbit=None):
        self._Orbit = Orbit

        # Data to query
        cycles = Orbit.Params.cycles
        self.states = self.update_dynamics(cycles)

        # Required to 'hold' satellite
        self._sample = None

    def query_states(self, tspan=(0,1)):
        ''' Convert to a common reference time and query state '''
        state_subspace = numpy.empty((6,0))
        
        step = self._Orbit.Params.step
        step_in_seconds = step*self._Orbit.Params.period_in_seconds()/self._Orbit.Params.period

        n_current = int(1+tspan[0]/step_in_seconds)
        n_future = int(1+tspan[-1]/step_in_seconds)
        nspan = (n_current, n_future)
        
        for n in nspan:        
            try:
                state = self._Orbit.Frames.barycenter_fixed_to_barycenter_lunar_inertial(self.states[:,n])
            except IndexError:
                total_samples_per_cycle = self._Orbit.Params.total_samples()
                n = n % total_samples_per_cycle
                state = self._Orbit.Frames.barycenter_fixed_to_barycenter_lunar_inertial(self.states[:,n])
            state_subspace = numpy.append(state_subspace, state, axis=1)
        return state_subspace

    def update_dynamics(self, cycles=1):
        ''' Update the vehicle dynamics '''
        return self._Orbit.get_barycenter_fixed_solution(cycles)
        