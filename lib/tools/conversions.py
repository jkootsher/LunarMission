import math
import numpy

from lib.tools.generators import xmat

def deg2rad(angle):
    ''' Degrees to radians '''
    return angle*math.pi/180

def rad2deg(angle):
    ''' Radians to degrees '''
    return angle*180/math.pi

def km2m(dist):
    ''' Kilometers to meters '''
    return dist*1000

def m2km(dist):
    ''' Meters to kilometers '''
    return dist/1000

def mrp2dcm(MRP):
    ''' MRP to DCM '''
    S_XMAT = xmat(MRP[:,-1])
    norm_squared = numpy.linalg.norm(MRP)**2
    R = numpy.identity(3) + (8*numpy.matmul(S_XMAT,S_XMAT) - 4*(1-norm_squared)*S_XMAT)/(1+norm_squared)**2
    return numpy.asarray(R)

def dcm2mrp(DCM):
    ''' DCM to MRP '''
    mrp = numpy.empty((3,1))
    if abs(numpy.trace(DCM)+1) < 1e-12:
        C = 1
    else:
        C = math.sqrt(numpy.trace(DCM)+1)
        C = 1/(C*(C+2))
    mrp[0,0] = C*(DCM[1,2]-DCM[2,1])
    mrp[1,0] = C*(DCM[2,0]-DCM[0,2])
    mrp[2,0] = C*(DCM[0,1]-DCM[1,0])
    return mrp

def unit_vector(v):
    ''' Convert to unit vector '''
    v_norm =numpy.linalg.norm(v)
    if v_norm < 1e-12:
        v = numpy.zeros(v.shape)
    else:
        v = v/v_norm
    return v