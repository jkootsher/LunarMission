import numpy
import matplotlib.pyplot as mplot

from lib.gnc.cr3bp import LCR3BP
from lib.gnc.kepler import Kepler
from lib.models.orbiter import Satellite
from lib.models.lunar_relay import Target
from lib.simcore.support.variables import GLOBALS
from lib.gnc.lunar.dynamics import Orbit as Lunar
from lib.gnc.libration.dynamics import Orbit as L2

ORBIT_CYCLES = 2

elements = Kepler(e=0.6, a=1000, i=56.2, argp=90, raan=0)
Lunar_Orbit = Lunar(elements)

elements = LCR3BP(h=2e-4, sbdy=GLOBALS.LUNAR['MU'], lbdy=GLOBALS.EARTH['MU'], cycles=ORBIT_CYCLES)
Halo_Orbit = L2(elements)
L2_Satellite = Target(Halo_Orbit)


t_min = 0
t_max = ORBIT_CYCLES*Lunar_Orbit.period

Orbiter = Satellite(Lunar_Orbit)
Orbiter.Target = L2_Satellite
vehicle_state = Orbiter.update_dynamics(delta=.1, dt=(t_min,t_max))

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