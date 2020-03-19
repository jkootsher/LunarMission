import math
import numpy

def gmat(r, mu):
    ''' Generate the matrix G that represents the orbital gravitation as a function of position '''

    r1_norm = math.sqrt((mu+r[0])**2 + r[1]**2 + r[2]**2)
    r2_norm = math.sqrt((r[0]-(1-mu))**2 + r[1]**2 + r[2]**2)

    u11 = 1 - (1-mu)*(1/r1_norm**3 - 3*(r[0]+mu)**2/r1_norm**5) - mu*(1/r2_norm**3 - 3*(r[0] - (1-mu))**2/r2_norm**5)
    u22 = 1 - (1-mu)*(1/r1_norm**3 - 3*(r[1]+mu)**2/r1_norm**5) - mu*(1/r2_norm**3 - 3*(r[1] - (1-mu))**2/r2_norm**5)
    u33 = -(1-mu)*(1/r1_norm**3 - 3*r[2]**2/r1_norm**5) - mu*(1/r2_norm**3 - 3*r[2]**2/r2_norm**5)

    u12 = 3*(1-mu)*r[1]*(r[0]+mu)/r1_norm**5 + 3*mu*r[1]*(r[0] - (1-mu))/r2_norm**5
    u13 = 3*(1-mu)*r[2]*(r[0]+mu)/r1_norm**5 + 3*mu*r[2]*(r[0] - (1-mu))/r2_norm**5

    u21 = u12
    u23 = 3*(1-mu)*r[1]*(r[2]+mu)/r1_norm**5 + 3*mu*r[1]*r[2]/r2_norm**5

    u31 = u13
    u32 = u23

    gmat = numpy.array([[u11,u12,u13], [u21, u22, u23], [u31, u32, u33]])
    return gmat
            