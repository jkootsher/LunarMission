import numpy
import matplotlib.pyplot as mplot

from lib.gnc.cr3bp import LCR3BP
from lib.gnc.kepler import Kepler
from lib.models.orbiter import Satellite
from lib.models.lunar_relay import Target
from lib.simcore.support.variables import GLOBALS
from lib.gnc.lunar.dynamics import Orbit as Lunar
from lib.gnc.libration.dynamics import Orbit as L2

ORBIT_CYCLES = 3

# Configure orbital parameters
elm_moon = Kepler(e=0.6, a=500, i=56.2, argp=90, raan=0)
elm_halo = LCR3BP(h=2e-4, sbdy=GLOBALS.LUNAR['MU'], lbdy=GLOBALS.EARTH['MU'], cycles=ORBIT_CYCLES)

# Required orbital physics
Halo_Orbit = L2(elm_halo)
Lunar_Orbit = Lunar(elm_moon)

# Create the satellites
Orbiter = Satellite(Lunar_Orbit)
L2_Satellite = Target(Halo_Orbit)

# Configure the orbiter tracking (optional)
Orbiter.configure_tracking(L2_Satellite)

# Calculate the orbiter dynamics
tspan = (0, ORBIT_CYCLES*Lunar_Orbit.period)
vehicle_state = Orbiter.update_dynamics(delta=.1, dt=tspan)


# Plots for analysis
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