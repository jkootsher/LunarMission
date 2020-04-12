import math
import numpy

import lib.tools.conversions as convert

from lib.gnc.kepler import Kepler
from lib.simcore.support.variables import GLOBALS
from lib.tools.math_toolbox.ode.solver import newton_raphson


class Orbit(object):
    ''' Dynamics for a Lunar Polar Orbit '''

    def __init__(self, Kepler=None):
        self._KOE = Kepler

        # Calculate Kepler orbit period
        self.MU = GLOBALS.LUNAR['MU']
        time_squared = Kepler.semi_major_axis**3/self.MU
        self.period = 2*math.pi*math.sqrt(time_squared)

    def get_pqw_solution(self, dt=0):
        ''' Compute the low lunar orbit (LLO) '''
        REVOLUTION = convert.deg2rad(GLOBALS.CONSTANTS['DEGREES_IN_CIRCLE'])
        
        # Guess at mean anomaly
        M = dt*(2*math.pi)/self.period

        # Find the eccentric anomaly
        eccentric_anomaly = newton_raphson(Kepler.mean_anomaly, df=Kepler.mean_anomaly_df, \
                                            x=M, e=self._KOE.eccentricity, M=M)

        # Normalize to a single rotation
        eccentric_anomaly = eccentric_anomaly - REVOLUTION*math.floor(eccentric_anomaly/REVOLUTION)

        # Calculate the true anomaly
        ratio = math.sqrt((1+self._KOE.eccentricity)/(1-self._KOE.eccentricity))
        self._KOE.true_anomaly = 2*math.atan(ratio*math.tan(eccentric_anomaly/2))

        # Perifocal state calculations
        ratio = (1-self._KOE.eccentricity**2)/(1+self._KOE.eccentricity*math.cos(self._KOE.true_anomaly))
        r_magnitude = self._KOE.semi_major_axis*ratio

        r_p = r_magnitude*math.cos(self._KOE.true_anomaly)
        r_q = r_magnitude*math.sin(self._KOE.true_anomaly)
        r_w = 0

        ratio = self.MU/(self._KOE.semi_major_axis*(1-self._KOE.eccentricity**2))
        v_p = -math.sqrt(ratio)*math.sin(self._KOE.true_anomaly)
        v_q = +math.sqrt(ratio)*(self._KOE.eccentricity + math.cos(self._KOE.true_anomaly))
        v_w = 0

        r_pqw = numpy.array([[r_p], [r_q], [r_w]])
        v_pqw = numpy.array([[v_p], [v_q], [v_w]])

        state_vector = numpy.append(r_pqw, v_pqw)
        state_vector = numpy.reshape(state_vector, (6,1))
        return state_vector
