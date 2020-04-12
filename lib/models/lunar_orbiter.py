from lib.gnc.lunar.frames import Frame
from lib.gnc.lunar.dynamics import Orbit

class Lun(object):
    ''' Lunar Model Class '''

    def __init__(self):
        self.state_vector = None

        self.frames = Frame()
