class GLOBALS(object):
    ''' Container for various globals and constants '''
    DEFAULT = {
        'MU' :    1,              # km^3/s^2
        'MASS':   1,              # kg
        'RADIUS': 1,              # km
    }

    LUNAR = {
        'MU':     4.90487e03,     # km^3/s^2
        'MASS':   0.07346e24,     # kg
        'RADIUS': 1738.10e00,     # km

        'TO_LIB2': 6.1347e03,     # km
        'T_ORBIT': 375.667e03,    # sec
    }

    EARTH = {
        'MU':     3.98600e05,     # km^3/s^2
        'MASS':   5.97240e24,     # kg
        'RADIUS': 6378.10e00,     # km

        'TO_MOON': 384.4e03,      # km
    }

    CONSTANTS = {
        'GRAV': 6.6743e-20,        # km^3/(kg*s^2)
        'VECTOR_SIZE': 3,          # n/a
        'STATE_VECTOR_DIM': 6,     # n/a
        'DEGREES_IN_CIRCLE': 360,  # n/a
    }

# Put the control torque and decay time elsewhere
    PARAMETERS = {
        'DECAY_TIME':  60,         # sec
        'CTRL_TORQUE': [0,0,0],    # Nm
        'SAMPLE_RATE': 1.0,        # Hz
    }
