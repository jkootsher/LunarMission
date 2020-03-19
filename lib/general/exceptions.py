class Verify(object):
    ''' Simple class for verifying data types '''
    
    @classmethod
    def is_numeric(cls, arg):
        ''' Check if arg is a valid number '''
        try:
            float(arg)
            return True
        except ValueError:
            return False


class ReferenceFrameError(Exception):
    def __init__(self):
        Exception.__init__(self, "Vectors are not in the same reference frame.")
        return

class VectorIndexError(Exception):
    def __init__(self):
        Exception.__init__(self, "Vectors are of indices [i, j, k] only.")
        return