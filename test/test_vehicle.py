from lib.simcore.support.variables import GLOBALS
from lib.models.lunar_orbiter import Satellite
from lib.gnc.libration.dynamics import Orbit as L2
from lib.gnc.lunar.dynamics import Orbit
from lib.gnc.kepler import Kepler

from lib.simcore.function_library import axisEqual3D

import math
import numpy
import matplotlib.pyplot as mplot

from mpl_toolkits.mplot3d import Axes3D

elements = Kepler(e=0.6, a=400, i=0, argp=0, raan=0)
Lunar_Orbit = Orbit(elements)

# t = 2e-4                   # step size
# Torbit = 3.1959             # multiply by 'T' to get time
# span_dt = (0,Torbit/2)      # time span for vehicle_state single orbit

# # TODO: CREATE A FUNCTION TO HANDLE DIFFERENT SAMPLING INTERVALS FOR COMM FRAME COMPARISON

# # HALO model (CRTBP)
# Halo_Orbit = L2(GLOBALS.LUNAR['MU'], GLOBALS.EARTH['MU'])

# # Change the default parameters (optional)
# position = (1.1503, 0, 0.1459)
# velocity = (0, -0.2180, 0)
# Halo_Orbit.change_initial_conditions(position, velocity)

# # Generate trajectory of HALO (N=1 default)
# trajectory = Halo_Orbit.get_barycenter_fixed_solution(N=4, delta=dt, tspan=span_dt)

# state_space_lci = numpy.empty((6,0))

# for n in range(Lunar_Orbit.total_orbital_samples()):
#     # Moon
#     state_vector_lci = Lunar_Orbit.inertial_solution(n)
#     state_space_lci = numpy.append(state_space_lci, state_vector_lci, axis=1)

t_min = 0
t_max = 4*Lunar_Orbit.period

Test_Sat = Satellite(Lunar_Orbit)
vehicle_state = Test_Sat.update_dynamics(delta=.1, dt=(t_min,t_max))

fig1 = mplot.figure()
ax = fig1.gca()
ax.plot(vehicle_state[0,:], color='blue', label='i-axis')
ax.plot(vehicle_state[1,:], color='red', label='j-axis')
ax.plot(vehicle_state[2,:], color='purple', label='k-axis')

mplot.title('Low Lunar Orbit (LCI)', fontsize=14)
ax.set_xlabel('Samples', fontsize=8)
ax.set_ylabel('MRP Set', fontsize=8)
ax.legend(loc='upper left')

fig2 = mplot.figure()
ax = fig2.gca()
ax.plot(vehicle_state[3,:], color='blue', label='i-axis')
ax.plot(vehicle_state[4,:], color='red', label='j-axis')
ax.plot(vehicle_state[5,:], color='purple', label='k-axis')

mplot.title('Low Lunar Orbit (LCI)', fontsize=14)
ax.set_xlabel('Samples', fontsize=8)
ax.set_ylabel('Body Rates', fontsize=8)
ax.legend(loc='upper left')

mplot.show()