from lib.gnc.kepler import Kepler
from lib.gnc.lunar.frames import Frame
from lib.gnc.lunar.dynamics import Orbit
from lib.simcore.support.variables import GLOBALS
from lib.simcore.function_library import axisEqual3D

import math
import numpy
import matplotlib.pyplot as mplot

from mpl_toolkits.mplot3d import Axes3D
# https://phys.org/news/2006-12-paradigm-lunar-orbits.html

elements = Kepler(e=0.6, a=6541.4, i=56.2, argp=90, raan=90)
lunar_orbit = Orbit(elements)
lunar_frame = Frame(elements)

state_space_lci = numpy.empty((6,0))
state_space_eci = numpy.empty((6,0))

MU = GLOBALS.LUNAR['MU']
T_LUNAR = 2*math.pi*math.sqrt(elements.semi_major_axis**3/MU)

orbit_period = int(T_LUNAR + 1)

for dt in range(orbit_period):
    state_vector_pqw = lunar_orbit.get_pqw_solution(dt=dt)

    # Moon
    state_vector_lci = lunar_frame.pqw2lci(state_vector_pqw)
    state_space_lci = numpy.append(state_space_lci, state_vector_lci, axis=1)
    
    # Earth
    state_vector_eci = lunar_frame.lci2eci(state_vector_lci)
    state_space_eci = numpy.append(state_space_eci, state_vector_eci, axis=1)

# Define Earth model
earth_radius = GLOBALS.EARTH['RADIUS']
distance_to_moon = GLOBALS.EARTH['TO_MOON']

# Define Moon model with the vehicle orbit
lunar_radius = GLOBALS.LUNAR['RADIUS']
foci_separation_magnitude = 2*elements.semi_major_axis*elements.eccentricity
x_second_focus = foci_separation_magnitude*numpy.cos(+elements.inclination)
z_second_focus = foci_separation_magnitude*numpy.sin(-elements.inclination)

# Sphere definition...
u = numpy.linspace(0, 2*math.pi, 100)
v = numpy.linspace(0, math.pi, 100)

# ...for Moon
x_lunar = lunar_radius*numpy.outer(numpy.cos(u), numpy.sin(v))
y_lunar = lunar_radius*numpy.outer(numpy.sin(u), numpy.sin(v))
z_lunar = lunar_radius*numpy.outer(numpy.ones(numpy.size(u)), numpy.cos(v))

# ...for Earth
x_earth = earth_radius*numpy.outer(numpy.cos(u), numpy.sin(v))
y_earth = earth_radius*numpy.outer(numpy.sin(u), numpy.sin(v))
z_earth = earth_radius*numpy.outer(numpy.ones(numpy.size(u)), numpy.cos(v))

# Earth Centered Inertial
fig1 = mplot.figure()
ax = fig1.gca(projection='3d')

# Earth
ax.plot_wireframe(x_earth, y_earth, z_earth, color='blue', rstride=8, cstride=8, label='Earth Surface')

# Orbit
ax.plot(state_space_eci[0,:], state_space_eci[1,:], state_space_eci[2,:], label='Lunar Orbit')

# Moon
ax.plot_wireframe(distance_to_moon+x_lunar, y_lunar, z_lunar, color='red', rstride=8, cstride=8, label='Lunar Surface')

axisEqual3D(ax)
mplot.title('Low Lunar Orbit (ECI)', fontsize=14)
ax.set_xlabel('ECI X (km)', fontsize=8)
ax.set_ylabel('ECI Y (km)', fontsize=8)
ax.set_zlabel('ECI Z (km)', fontsize=8)
ax.legend(loc='upper left')


# Lunar Centered Inertial
fig2 = mplot.figure()
ax = fig2.gca(projection='3d')

# Orbit
ax.plot(state_space_lci[0,:], state_space_lci[1,:], state_space_lci[2,:], label='Lunar Orbit')

# Moon
ax.scatter(0, 0, 0, marker='o', color='green', label='Lunar CM (Interior Focus)')
ax.scatter(x_second_focus, 0, z_second_focus, marker='o', color='purple', label='Exterior Focus')
ax.plot_wireframe(x_lunar, y_lunar, z_lunar, color='red', rstride=8, cstride=8, label='Lunar Surface')

axisEqual3D(ax)
mplot.title('Low Lunar Orbit (LCI)', fontsize=14)
ax.set_xlabel('LCI X (km)', fontsize=8)
ax.set_ylabel('LCI Y (km)', fontsize=8)
ax.set_zlabel('LCI Z (km)', fontsize=8)
ax.legend(loc='upper left')
mplot.show()