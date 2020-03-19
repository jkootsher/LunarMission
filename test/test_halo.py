import numpy
import matplotlib.pyplot as mplot
from lib.gnc.orbits import dynamics
from mpl_toolkits.mplot3d import Axes3D
from lib.general.variables import GLOBALS

# T = GLOBALS.LUNAR['ORBIT_T'] conversion factor

# SUCCESSFUL TEST RUN UNDER THE FOLLOWING CONDITIONS:
#
# dynamics.location = (1.0874, 0, .2022)
# dynamics.velocity = (0, -.1926, 0)
# Torbit = 2.4499
# dt = 2e-5           order of 1e-5 gives closed path
#
# The above parameters for location and velocity
# are the defaults in the gnc.orbits.dynamics.Halo class
#
# Initial conditions taken from Table 3.5 on page 71
# https://engineering.purdue.edu/people/kathleen.howell.1/Publications/Masters/2006_Grebow.pdf
# L2 Lyapunov Orbital Parameters for Northern L2 HALO family
#

# Timekeeping
dt = 2e-4                   # step size
Torbit = 3.1959             # multiply by 'T' to get time
tfinal = Torbit/2           # half the orbital period
samples = int(tfinal/dt)    # number of samples
span_dt = numpy.empty((0))  # initialize the time range

for n in range(samples):
    span_dt = numpy.append(span_dt, [n])

# HALO model (CRTBP)
halo = dynamics.Halo(GLOBALS.LUNAR['MU'], GLOBALS.EARTH['MU'])

# Change the default parameters (optional)
position = (1.1503, 0, 0.1459)
velocity = (0, -0.2180, 0)
halo.change_initial_conditions(position, velocity)

# Locations of interest
L2_xlocation = halo.L2
Earth_xlocation = -halo.mu
Lunar_xlocation = 1-halo.mu

# Generate trajectory of HALO (N=1 Lyapunov)
trajectory = halo.compute(step_size=dt, samples=span_dt)

# Plot HALO model
fig = mplot.figure()
ax = fig.gca(projection='3d')
ax.plot(trajectory[0,:], trajectory[1,:], trajectory[2,:])
ax.scatter(L2_xlocation, 0, 0, marker='o', color='green')
ax.scatter(Earth_xlocation, 0, 0, marker='v', color='purple')
ax.scatter(Lunar_xlocation, 0, 0, marker='^', color='red')
mplot.show()