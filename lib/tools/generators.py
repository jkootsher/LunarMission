import math
import numpy


def gmat(r, mu):
    ''' Generate the matrix G that represents the orbital gravitation as a function of position '''

    r1_norm = math.sqrt((r[0]+mu)**2 + r[1]**2 + r[2]**2)
    r2_norm = math.sqrt((r[0]-(1-mu))**2 + r[1]**2 + r[2]**2)

    u11 = 1 - (1-mu)*(1/r1_norm**3 - 3*(r[0]+mu)**2/r1_norm**5) - mu*(1/r2_norm**3 - 3*(r[0] - (1-mu))**2/r2_norm**5)
    u22 = 1 - (1-mu)*(1/r1_norm**3 - 3*r[1]**2/r1_norm**5) - mu*(1/r2_norm**3 - 3*r[1]**2/r2_norm**5)
    u33 = 0 - (1-mu)*(1/r1_norm**3 - 3*r[2]**2/r1_norm**5) - mu*(1/r2_norm**3 - 3*r[2]**2/r2_norm**5)

    u12 = 3*(1-mu)*r[1]*(r[0]+mu)/r1_norm**5 + 3*mu*r[1]*(r[0] - (1-mu))/r2_norm**5
    u13 = 3*(1-mu)*r[2]*(r[0]+mu)/r1_norm**5 + 3*mu*r[2]*(r[0] - (1-mu))/r2_norm**5

    u21 = u12
    u23 = 3*(1-mu)*r[1]*r[2]/r1_norm**5 + 3*mu*r[1]*r[2]/r2_norm**5

    u31 = u13
    u32 = u23

    gmat = numpy.array([[u11,u12,u13], [u21,u22,u23], [u31,u32,u33]])
    return gmat

def xmat(s):
    ''' Generate the matrix X that represents a cross product from s '''
    xmat = numpy.array([[0, -s[2], s[1]], [s[2], 0, -s[0]], [-s[1], s[0], 0]])
    return xmat

def euler_313(phi):
    ''' Generate the 313 rotation sequence (angles must be radian) '''
    r11 = math.cos(phi[2])*math.cos(phi[0]) - math.sin(phi[2])*math.cos(phi[1])*math.sin(phi[0])
    r12 = -math.cos(phi[2])*math.sin(phi[0]) - math.sin(phi[2])*math.cos(phi[1])*math.cos(phi[0])
    r13 = math.sin(phi[2])*math.sin(phi[1])
    r21 = math.sin(phi[2])*math.cos(phi[0]) + math.cos(phi[2])*math.cos(phi[1])*math.sin(phi[0])
    r22 = -math.sin(phi[2])*math.sin(phi[0]) + math.cos(phi[2])*math.cos(phi[1])*math.cos(phi[0])
    r23 = -math.cos(phi[2])*math.sin(phi[1])
    r31 = math.sin(phi[1])*math.sin(phi[0])
    r32 = math.sin(phi[1])*math.cos(phi[0])
    r33 = math.cos(phi[1])
    return numpy.asarray([[r11, r12, r13], [r21, r22, r23], [r31, r32, r33]])

def euler_232(phi):
    ''' Generate the 232 rotation sequence (angles must be radian) '''
    r11 = math.cos(phi[2])*math.cos(phi[1])*math.cos(phi[0]) - math.sin(phi[2])*math.sin(phi[0])
    r12 = math.cos(phi[2])*math.sin(phi[1])
    r13 = -math.sin(phi[2])*math.cos(phi[0]) - math.cos(phi[2])*math.cos(phi[1])*math.sin(phi[0])
    r21 = -math.sin(phi[1])*math.cos(phi[0])
    r22 = math.cos(phi[1])
    r23 = math.sin(phi[1])*math.sin(phi[0])
    r31 = math.cos(phi[2])*math.sin(phi[0]) + math.sin(phi[2])*math.cos(phi[1])*math.cos(phi[0])
    r32 = math.sin(phi[2])*math.sin(phi[1])
    r33 = math.cos(phi[2])*math.cos(phi[0]) - math.sin(phi[2])*math.cos(phi[1])*math.sin(phi[0])
    return numpy.asarray([[r11, r12, r13], [r21, r22, r23], [r31, r32, r33]])