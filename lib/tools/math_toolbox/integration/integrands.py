import math
import numpy

from lib.tools.generators import gmat
from lib.tools.generators import xmat
from lib.simcore.support.variables import GLOBALS

class CRTBP(object):

    @classmethod
    def vfhandle(cls, **kwargs):
        ''' Function handle (state vector representation) '''
        mu = kwargs['mu']

        state = numpy.empty((0))
        state = numpy.append(state, kwargs['y'])

        r1_norm = math.sqrt((mu+state[0])**2 + state[1]**2 + state[2]**2)
        r2_norm = math.sqrt((1-mu-state[0])**2 + state[1]**2 + state[2]**2)

        mass1 = 1 - mu
        mass2 = mu

        G = 1

        df1 = state[3]
        df2 = state[4]
        df3 = state[5]
        df4 = state[0] + 2*state[4] - G*mass1*(mu+state[0])/r1_norm**3 + G*mass2*(1-mu-state[0])/r2_norm**3
        df5 = state[1] - 2*state[3] - G*mass1*state[1]/r1_norm**3 - G*mass2*state[1]/r2_norm**3
        df6 = -G*mass1*state[2]/r1_norm**3 - G*mass2*state[2]/r2_norm**3
        df = numpy.asarray([[df1, df2, df3, df4, df5, df6]])

        return df.T

    @classmethod
    def mfhandle(cls, **kwargs):
        ''' Function handle (state matrix representation) '''
        # This is the ODE A' = dU*A
        mu = kwargs['mu']

        state = numpy.empty((0))
        state = numpy.append(state, kwargs['y'])
 
        dim_3 = GLOBALS.CONSTANTS['VECTOR_SIZE']
        dim_6 = GLOBALS.CONSTANTS['STATE_VECTOR_DIM']
        
        # Create partitions
        O = numpy.zeros((dim_3,dim_3))
        I = numpy.identity(dim_3)
        G = gmat(state[36:39], mu)
        K = numpy.zeros((dim_3,dim_3))

        K[0,1] = 2
        K[1,0] = -K[0,1]

        # Each partition is 3x3 so matrix is 6x6
        # dU = | O  I |
        #      | G  K |
        dU = numpy.zeros((dim_6,dim_6))
        dU[:,:dim_3] = numpy.concatenate((O, G), axis=0)
        dU[:,dim_3:] = numpy.concatenate((I, K), axis=0)

        # Build the A matrix from state data (6x6)
        A = numpy.reshape(state[:dim_6**2], (dim_6,dim_6))

        # ODE A' = dU*A
        dA = numpy.matmul(dU, A)
        dA = numpy.reshape(dA, (dim_6**2,1))        

        # ODE x' = f(x)
        kwargs['y'] = state[36:]
        df = numpy.asarray(cls.vfhandle(**kwargs))
        dF = numpy.append(dA, df, axis=0)

        return dF


class MRP(object):

    @classmethod
    def fhandle(cls, **kwargs):
        ''' Function handle for equation of motion '''
        state = kwargs['y']
        moments = kwargs['moments']
        control_torque = kwargs['control_torque']

        # Break up state vector into MRPs and body rates
        s = state[0:3]
        w = state[3:6]

        # Cross product matrix representation
        S_XMAT = xmat(s[:,-1])
        W_XMAT = xmat(w[:,-1])

        # Equation of motion
        S_MAT = .25*((1-numpy.linalg.norm(s)**2)*numpy.identity(3) + 2*S_XMAT + 2*numpy.dot(s,s.T))
        s = numpy.matmul(S_MAT, w)

        # Rate calculation
        I_INV = numpy.linalg.inv(moments)
        W_MAT = numpy.matmul(W_XMAT, moments)
        W_MAT = numpy.matmul(I_INV, W_MAT)
        w = -numpy.matmul(W_MAT, w) + control_torque

        # Updated state
        state = numpy.append(s, w, axis=0)
        return state
