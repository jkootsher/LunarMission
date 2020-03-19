import math

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