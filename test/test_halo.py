import numpy
import matplotlib.pyplot as mplot

from mpl_toolkits.mplot3d import Axes3D

from lib.gnc.cr3bp import LCR3BP
from lib.gnc.libration.frames import Frames
from lib.gnc.libration.dynamics import Orbit
from lib.simcore.support.variables import GLOBALS

# T = GLOBALS.LUNAR['ORBIT_T'] conversion factor

# SUCCESSFUL TEST RUN UNDER THE FOLLOWING CONDITIONS:
#
# dynamics.location = (1.0874, 0, .2022)
# dynamics.velocity = (0, -.1926, 0)
# Torbit = 2.4499
# dt = 2e-4           order of 1e-4 gives closed path
#
# The above parameters for location and velocity
# are the defaults in the gnc.orbits.dynamics.Halo class
#
# Initial conditions taken from Table 3.5 on page 58 (pdf 71)
# https://engineering.purdue.edu/people/kathleen.howell.1/Publications/Masters/2006_Grebow.pdf
# L2 Lyapunov Orbital Parameters for Northern L2 HALO family
#

# HALO model (CRTBP)
elements = LCR3BP(h=2e-4, sbdy=GLOBALS.LUNAR['MU'], lbdy=GLOBALS.EARTH['MU'])
Halo_Orbit = Orbit(elements)

# Generate trajectory of HALO (N=1 default)
trajectory = Halo_Orbit.get_barycenter_fixed_solution(N=4)

# Locations of interest
L2_xlocation = Halo_Orbit.L2
Earth_xlocation = -Halo_Orbit.mu
Lunar_xlocation = 1-Halo_Orbit.mu

# ECI frame
trajectory_eci = numpy.empty((6,0))

elapsed_time = trajectory.shape[1]
for t in range(elapsed_time):
    eci_vector = Halo_Orbit.Frames.barycenter_fixed_to_barycenter_earth_inertial(trajectory[:,t])
    trajectory_eci = numpy.append(trajectory_eci, eci_vector, axis=1)
    
time_in_eci_frame = [idx*elements.step*GLOBALS.LUNAR['T_ORBIT'] for idx in range(trajectory_eci.shape[1])]

# Verify the satellite velocity normal to the orbital plane is
# approx. the mean orbital speed of moon (1.022 km/s). This should
# occur when the satellite is directly above or below the moon (~1{3}/4 T)
lunar_match_time = int(elements.period/elements.step)//5
print(trajectory_eci[4,lunar_match_time])

# Plot HALO model
fig1 = mplot.figure()
ax = fig1.gca(projection='3d')
ax.plot(trajectory[0,:], trajectory[1,:], trajectory[2,:])
ax.scatter(L2_xlocation, 0, 0, marker='o', color='green', label='L2')
ax.scatter(Earth_xlocation, 0, 0, marker='v', color='purple', label='Earth')
ax.scatter(Lunar_xlocation, 0, 0, marker='^', color='red', label='Moon')
mplot.title('L2 Halo Trajectory', fontsize=14)
ax.set_xlabel('3-Body Barycentric Frame (i-axis)', fontsize=8)
ax.set_ylabel('3-Body Barycentric Frame (j-axis)', fontsize=8)
ax.set_zlabel('3-Body Barycentric Frame (k-axis)', fontsize=8)
ax.legend(loc='upper left')

fig2 = mplot.figure()
ax = fig2.gca()
ax.plot(trajectory[0,:], color='blue', label='i-axis')
ax.plot(trajectory[1,:], color='red', label='j-axis')
ax.plot(trajectory[2,:], color='purple', label='k-axis')
mplot.title('L2 Halo Trajectory (Barycentric Frame)', fontsize=14)
ax.set_xlabel('Time (common barycentric unit)', fontsize=8)
ax.set_ylabel('Distance (common barycentric units)', fontsize=8)
ax.legend(loc='upper left')

fig3 = mplot.figure()
ax = fig3.gca()
barycentric_common_time = ax.get_xticks()
eci_frame_ref_time = barycentric_common_time*GLOBALS.LUNAR['T_ORBIT']*elements.step*trajectory_eci.shape[1]
ax.set_xticks(eci_frame_ref_time)
ax.plot(time_in_eci_frame, trajectory_eci[0,:], color='blue', label='x-axis')
ax.plot(time_in_eci_frame, trajectory_eci[1,:], color='red', label='y-axis')
ax.plot(time_in_eci_frame, trajectory_eci[2,:], color='purple', label='z-axis')
mplot.title('L2 Halo Trajectory (ECI Frame)', fontsize=14)
ax.set_xlabel('Time (sec)', fontsize=8)
ax.set_ylabel('Distance (km)', fontsize=8)
ax.legend(loc='upper left')

fig4 = mplot.figure()
ax = fig4.gca()
ax.plot(trajectory[3,:], color='blue', label='i-axis')
ax.plot(trajectory[4,:], color='red', label='j-axis')
ax.plot(trajectory[5,:], color='purple', label='k-axis')
mplot.title('L2 Halo Velocities (Barycentric Frame)', fontsize=14)
ax.set_xlabel('Time (common barycentric unit)', fontsize=8)
ax.set_ylabel('Velocity (common barycentric units)', fontsize=8)
ax.legend(loc='upper left')

fig5 = mplot.figure()
ax = fig5.gca()
barycentric_common_time = ax.get_xticks()
eci_frame_ref_time = barycentric_common_time*GLOBALS.LUNAR['T_ORBIT']*elements.step*trajectory_eci.shape[1]
ax.set_xticks(eci_frame_ref_time)
ax.plot(time_in_eci_frame, trajectory_eci[3,:], color='blue', label='x-axis')
ax.plot(time_in_eci_frame, trajectory_eci[4,:], color='red', label='y-axis')
ax.plot(time_in_eci_frame, trajectory_eci[5,:], color='purple', label='z-axis')
mplot.title('L2 Halo Velocities (ECI Frame)', fontsize=14)
ax.set_xlabel('Time (sec)', fontsize=8)
ax.set_ylabel('Velocity (km/s)', fontsize=8)
ax.legend(loc='upper left')
mplot.show()