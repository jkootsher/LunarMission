import numpy

def feval(fhandle=None, **kwargs):
    ''' Evaluate function fhandle at args '''
    return fhandle(**kwargs)


def axisEqual3D(ax):
    ''' Auto correct the 3D axis aspect ratio '''
    # Source: https://stackoverflow.com/questions/8130823/set-matplotlib-3d-plot-aspect-ratio
    extents = numpy.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    sz = extents[:,1] - extents[:,0]
    centers = numpy.mean(extents, axis=1)
    maxsize = max(abs(sz))
    r = maxsize/2
    for ctr, dim in zip(centers, 'xyz'):
        getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)
    return