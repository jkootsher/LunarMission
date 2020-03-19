import numpy
from lib.general.variables import GLOBALS
from lib.gnc.orbits import dynamics
import matplotlib.pyplot as mplot
import lib.tools.math_toolbox.ode.solver as solver

h=.01
span_dt = numpy.empty((0))
for n in range(1000):
    span_dt = numpy.append(span_dt, [n])


S = solver.test(step_size=h, state=[3], samples=span_dt)
mplot.plot(S)
mplot.show()