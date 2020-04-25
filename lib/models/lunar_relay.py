import numpy


class Target(object):
    ''' Generic Relay Target Orbiter '''

    def __init__(self, Orbit=None):
        self._Orbit = Orbit

        # Data to query
        cycles = Orbit.Params.cycles
        self.states = self.update_dynamics(cycles)

    def query_states(self, tspan=(0,1)):
        ''' Convert to a common reference time and query state '''
        state_subspace = numpy.empty((6,0))
        
        step = self._Orbit.Params.step
        nspan = (int(1+tspan[0]/step), int(1+tspan[-1]/step))
        
        for n in nspan:
            try:
                state = self._Orbit.Frames.barycenter_fixed_to_barycenter_lunar_inertial(self.states[:,n])
            except IndexError:
                total_samples_per_cycle = self._Orbit.Params.total_samples()
                n = total_samples_per_cycle % n
                state = self._Orbit.Frames.barycenter_fixed_to_barycenter_lunar_inertial(self.states[:,n])
            state_subspace = numpy.append(state_subspace, state, axis=1)
        return state_subspace

    def update_dynamics(self, cycles=1):
        ''' Update the vehicle dynamics '''
        return self._Orbit.get_barycenter_fixed_solution(cycles)
        