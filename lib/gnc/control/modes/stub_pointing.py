import numpy


class Stub(object):
    ''' Controller Stub '''

    def relative_rates(self, **kwargs):
        ''' Passthrough '''
        return numpy.ones((6,1))

    def inertial_to_control(self, **kwargs):
        ''' Passthrough '''
        return numpy.identity(3)