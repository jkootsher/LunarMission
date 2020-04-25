import numpy


class Communication(Pointing):
    ''' Communications Pointing Control Mode '''

    def __init__(self, Orbit=None):
        self.Orbit = Orbit

    def relative_rates(self, delta_r_now=None, delta_r_next=None):
        ''' Communication frame dynamics between the lunar orbiter and the L2 satellite '''
        BYI2COM_now = self.inertial_to_control(delta_r_now)
        BYI2COM_next = self.inertial_to_control(delta_r_next)

        BYI2COM_DT = BYI2COM_next - BYI2COM_now
        RATE_XMAT = -numpy.matmul(BYI2COM_now.T, BYI2COM_DT)

        rates = numpy.array([[RATE_XMAT[2,1], RATE_XMAT[0,2], RATE_XMAT[1,0]]])
        state_vector = numpy.append(delta_r_now, rates, axis=0)
        return state_vector

    def inertial_to_control(self, state_vector=None):
        ''' Inertial origin frame to orbiter communication frame '''
        displacement = state_vector[0:3]
        norm = numpy.linalg.norm(displacement)

        BYI2COM = numpy.zeros((3,3))
        BYI2COM[0,:] = -displacement/norm

        n_k = numpy.array([[0, 0, 1]])
        row_2 = numpy.cross(displacement, n_k)
        row_2 = row_2/numpy.linalg.norm(row_2)
        BYI2COM[1,:] = row_2
        BYI2COM[2,:] = numpy.cross(BYI2COM[0,:], BYI2COM[1,:])
        return BYI2COM

    def control_to_inertial(self, state_vector=None):
        ''' Orbiter communication frame to inertial origin frame '''
        BYI2COM = self.inertial_to_control(state_vector)
        return BYI2COM.T