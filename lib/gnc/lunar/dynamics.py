import math
import numpy

import lib.tools.conversions as convert

from lib.gnc.lunar.frames import Frames
from lib.tools.math_toolbox.ode import solver
from lib.simcore.support.variables import GLOBALS


class Orbit(object):
    ''' Dynamics for a Lunar Polar Orbit '''

    def __init__(self, Kepler=None, phase=0):
        self.Kepler = Kepler
        self.Frames = Frames(Kepler)

        # Set the orbital phase
        self.phase = phase

        # Calculate Kepler orbit period
        self.mu = GLOBALS.LUNAR['MU']
        time_squared = Kepler.semi_major_axis**3/self.mu
        self.period = 2*math.pi*math.sqrt(time_squared)

    def total_orbital_samples(self):
        ''' Get the total number of orbit samples from the period '''
        return int(self.period+1)

    def pqw_solution(self, n=0):
        ''' Compute the specified lunar orbit from Kepler parameters starting at the given phase '''
        revolution = convert.deg2rad(GLOBALS.CONSTANTS['DEGREES_IN_CIRCLE'])
        phase = convert.deg2rad(self.phase)

        # Guess at mean anomaly
        M = (n+phase)*(2*math.pi)/self.period

        # Find the eccentric anomaly
        eccentric_anomaly = solver.newton_raphson(self.Kepler.mean_anomaly,         \
                                                df=self.Kepler.mean_anomaly_df,     \
                                                x=M, e=self.Kepler.eccentricity, M=M)

        # Normalize to a single rotation
        eccentric_anomaly = eccentric_anomaly - revolution*math.floor(eccentric_anomaly/revolution)

        # Calculate the true anomaly
        ratio = math.sqrt((1+self.Kepler.eccentricity)/(1-self.Kepler.eccentricity))
        self.Kepler.true_anomaly = 2*math.atan(ratio*math.tan(eccentric_anomaly/2))

        # Perifocal state calculations
        ratio = (1-self.Kepler.eccentricity**2)/(1+self.Kepler.eccentricity*math.cos(self.Kepler.true_anomaly))
        r_magnitude = self.Kepler.semi_major_axis*ratio

        r_p = r_magnitude*math.cos(self.Kepler.true_anomaly)
        r_q = r_magnitude*math.sin(self.Kepler.true_anomaly)
        r_w = 0

        ratio = self.mu/(self.Kepler.semi_major_axis*(1-self.Kepler.eccentricity**2))
        v_p = -math.sqrt(ratio)*math.sin(self.Kepler.true_anomaly)
        v_q = +math.sqrt(ratio)*(self.Kepler.eccentricity + math.cos(self.Kepler.true_anomaly))
        v_w = 0

        r_pqw = numpy.array([[r_p], [r_q], [r_w]])
        v_pqw = numpy.array([[v_p], [v_q], [v_w]])

        # Generate the solution set
        state_vector = numpy.append(r_pqw, v_pqw)
        state_vector = numpy.reshape(state_vector, (6,1))
        return state_vector

    def inertial_solution(self, n=0):
        ''' Compute the inertial solution for the specified orbit '''
        state_vector = self.pqw_solution(n)
        PQW2BYI = self.Frames.perifocal_to_barycenter_lunar_inertial()
        return numpy.matmul(PQW2BYI, state_vector)
