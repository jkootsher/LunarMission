import numpy
import lib.tools.math_toolbox.ode.solver as solver

from lib.general.datatypes import Vector
from lib.general.variables import GLOBALS
from lib.tools.physics_toolbox import energy
from lib.tools.math_toolbox.integration.integrands import CRTBP as Functional


class Halo(object):
    ''' Dynamics for the L2 Halo Orbit '''
    def __init__(self, small_bdy=GLOBALS.DEFAULT['MU'], large_bdy=GLOBALS.DEFAULT['MU']):
        self.mu = small_bdy/(small_bdy + large_bdy)
        self.position = Vector((1.0874,0.0,0.2020), 'L2 Inertial')
        self.velocity = Vector((0.0,-0.2054,0.0), 'L2 Inertial')

        self.L2 = self._lagrange_point()
        return

    def _lagrange_point(self, true_root=0):
        ''' Calculate the L2 libration (langrange) point location '''
        mu_a = self.mu
        mu_b = 1 - mu_a

        # Taylor expansion
        poly_coeffs = [1.0,                                   \
                       2*(mu_a-mu_b),                         \
                       mu_b**2 - 4*mu_a*mu_b + mu_a**2,       \
                       2*mu_a*mu_b*(mu_b-mu_a) - (mu_b+mu_a), \
                       (mu_a*mu_b)**2 + 2*(mu_b**2-mu_a**2),  \
                       -(mu_a**3+mu_b**3)]

        # Obtain the minimum 'largest' on the line formed
        # by the Earth-Moon system. This is the L2 point
        roots = numpy.roots(poly_coeffs)
        for idx in range(len(poly_coeffs)-1):
            if roots[idx] > -mu_a and roots[idx] > mu_b:
                true_root = numpy.absolute(roots[idx])
        return true_root

    def change_initial_conditions(self, position=None, velocity=None):
        ''' Change the initial orbital parameters '''
        self.position = Vector(position, 'L2 Inertial')
        self.velocity = Vector(velocity, 'L2 Inertial')
        return

    def compute(self, N=1, **kwargs):
        ''' Compute the orbit using N patched trajectories '''
        h = kwargs['step_size']
        state_space = numpy.empty((0))
        dimension = GLOBALS.CONSTANTS['STATE_VECTOR_DIM']
        current_state = list(self.position.at(0)) + list(self.velocity.at(0))
        current_state = numpy.asarray([current_state])

        for idx in range(N):
            kwargs['mu'] = self.mu
            kwargs['state'] = current_state.T
            [samples, vstate] = solver.ode113v(**kwargs)
            
            # Using the current state solution, calculate the
            # state matrix required for the current trajectory
            multibody_state = numpy.identity(dimension)
            multibody_state = numpy.reshape(multibody_state, (dimension**2,1))
            kwargs['state'] = numpy.append(multibody_state, current_state.T, axis=0)
            [samples, mstate]  = solver.ode113m(**kwargs)

            # State matrix vector representation conversion back to matrix
            A = numpy.reshape(mstate[:36,-1], (dimension,dimension))

            # Obtain the directional vector field
            kwargs['y'] = vstate
            vfield = Functional.vfhandle(**kwargs)

            # Solution to state matrix confined to vector field
            # generated by the corresponding state solution at time n
            dF = numpy.array([[A[3,0], A[3,4], vfield[3]], \
                              [A[5,0], A[5,4], vfield[5]], \
                              [A[1,0], A[1,4], vfield[1]]])

            inv_dF = numpy.linalg.inv(dF)
            push_state = numpy.asarray([[vstate[3,-1], vstate[5,-1], vstate[1,-1]]])
            
            # Estimate the initial conditions for the next trajectory
            # patch given the current directional state space solution
            next_state_estimate = numpy.asarray([[current_state[0,0], current_state[0,4], len(samples)*h]])
            next_state_estimate = next_state_estimate.T - numpy.matmul(inv_dF, push_state.T)            
            next_state = numpy.array([[next_state_estimate[0,0], 0, self.position.at(0)[2], 0, next_state_estimate[1,0], 0]])

            # Calculate an initial guess at the orbital
            # period for the next trajectory patch N
            time_estimate = int(next_state_estimate[2]/h)
            samples = [n for n in range(time_estimate)]
            kwargs['samples'] = numpy.asarray(samples)
            current_state = next_state
        
        # Orbital symmetry allows multiplication of the total time
        # by two so as to find the complete orbital path (initial time was T/2)
        samples = [n for n in range(2*len(kwargs['samples']))]
        kwargs['state'] = current_state.T
        kwargs['samples'] = numpy.asarray(samples)
        [deltaT,trajectory] = solver.ode113v(**kwargs)
        return trajectory