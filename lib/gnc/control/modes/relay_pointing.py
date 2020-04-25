import math
import numpy

from lib.gnc.control.frames import Frames
from lib.tools.conversions import unit_vector


class Communication(Frames):
    ''' Communications Pointing Control Mode '''

    def __init__(self):
        super(Communication, self).__init__()
        self.fname = 'communication'

    def track(self, Target=None, **kwargs):
        ''' Sync to the vehicle that is being tracked '''
        (ni,nf) = kwargs['nspan']
        tspan = (ni*kwargs['delta'],nf*kwargs['delta'])
        tracked_positions = Target.query_states(tspan)
        current_positions = kwargs['local']

        delta_positions = tracked_positions-current_positions
        return delta_positions

    def update_span(self, **kwargs):
        ''' Update the vehicle delta span information '''
        delta_position = kwargs['target'][0:3,0]
        vehicle_position = kwargs['local'][0:3,0]

        target_position = delta_position + vehicle_position

        position_norm = numpy.linalg.norm(vehicle_position)
        tracking_norm = numpy.linalg.norm(target_position)

        # Find conic section angle
        projection = numpy.dot(vehicle_position, target_position)
        delta_angle = math.acos(projection/(position_norm*tracking_norm))
        return delta_angle

    def relative_rates(self, **kwargs):
        ''' Relative rate calculation for the relay frame '''
        delta_present = kwargs['target'][:,0]
        delta_tracked = kwargs['target'][:,-1]        

        # Prepare the current state attitude
        kwargs['target'] = numpy.reshape(delta_present, (6,1))
        BYI2COM_n = self.inertial_to_control(**kwargs)

        # Prepare the future state attitude
        kwargs['target'] = numpy.reshape(delta_tracked, (6,1))
        BYI2COM_f = self.inertial_to_control(**kwargs)

        # Rate matrix (from derivative)
        DELTA = BYI2COM_f-BYI2COM_n
        W_XMAT = -numpy.matmul(BYI2COM_n.T, DELTA)

        w0 = W_XMAT[2,1]
        w1 = W_XMAT[0,2]
        w2 = W_XMAT[1,0]
        rates = numpy.asarray([w0,w1,w2])
        return numpy.reshape(rates, (3,1))

    def inertial_to_control(self, **kwargs):
        ''' Inertial origin frame to orbiter communication frame '''
        position = kwargs['target'][0:3,0]

        # Create axis along node line
        n = numpy.asarray([0,0,1])
        n_position = numpy.cross(position.T, n.T)

        # Construct the frame
        BYI2COM = numpy.zeros((3,3))
        BYI2COM[0,:] = -unit_vector(position.T)
        BYI2COM[1,:] = +unit_vector(n_position.T)
        BYI2COM[2,:] = numpy.cross(BYI2COM[0,:], BYI2COM[1,:])
        return BYI2COM

    def control_to_inertial(self, **kwargs):
        ''' Orbiter communication frame to inertial origin frame '''
        BYI2COM = self.inertial_to_control(**kwargs)
        return BYI2COM.T