import numpy
from lib.simcore.function_library import feval


def RK4(fhandle=None, **kwargs):
    ''' Runge-Kutta (4th order) numerical integrator '''
    h = kwargs['h']
    
    x_span = kwargs['x']
    x_delta = int((x_span[-1] - x_span[0])/h)
    x_initial = x_span[0]

    y_initial = kwargs['y']
    state_size = len(y_initial)

    x_update = x_initial
    y_update = y_initial

    x_state = numpy.asarray([x_update])
    y_state = numpy.asarray(y_update)
    
    for idx in range(x_delta):
        kwargs['x'] = x_update
        kwargs['y'] = y_update
        k1 = feval(fhandle, **kwargs)

        kwargs['x'] = x_update + h/2
        kwargs['y'] = y_update + k1*(h/2)
        k2 = feval(fhandle, **kwargs)

        kwargs['x'] = x_update + h/2
        kwargs['y'] = y_update + k2*(h/2)
        k3 = feval(fhandle, **kwargs)

        kwargs['x'] = x_update + h
        kwargs['y'] = y_update + k3*h
        k4 = feval(fhandle, **kwargs)

        x_update = x_update + h
        x_state = numpy.append(x_state, x_update)

        y_update = y_update + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
        y_state = numpy.append(y_state, y_update, axis=1)
        
    return [x_state, y_state]


def ABM4(fhandle=None, **kwargs):
    ''' Adam's & Bashforth method (4th order) for numerical integration '''
    h = kwargs['h']
    
    x_span = kwargs['x']
    x_delta = int((x_span[-1] - x_span[0])/h)
    x_initial = [x_span[0] + n*h for n in range(x_delta+1)]

    y_initial = kwargs['y']
    state_size = len(y_initial)

    kwargs['x'] = (x_initial[0], x_initial[3])
    kwargs['y'] = y_initial

    [x_update, y_update] = RK4(fhandle, **kwargs)
    
    x_state = numpy.asarray(x_update)
    y_state = numpy.asarray(y_update)

    for idx in range(3, x_delta):
        kwargs['x'] = x_update[idx]
        kwargs['y'] = y_update[:,idx]
        df0 = feval(fhandle, **kwargs)

        kwargs['x'] = x_update[idx-1]
        kwargs['y'] = y_update[:,idx-1]
        df1 = feval(fhandle, **kwargs)

        kwargs['x'] = x_update[idx-2]
        kwargs['y'] = y_update[:,idx-2]
        df2 = feval(fhandle, **kwargs)

        kwargs['x'] = x_update[idx-3]
        kwargs['y'] = y_update[:,idx-3]
        df3 = feval(fhandle, **kwargs)

        kwargs['x'] = x_update[idx] + h
        kwargs['y'] = y_update[:,idx] + (h/24)*(55*df0 - 59*df1 + 37*df2 - 9*df3)
        df = feval(fhandle, **kwargs)

        x_update = x_update[idx] + h
        x_state = numpy.append(x_state, x_update)

        y_update = y_update[:,[idx]] + (h/24)*(9*df + 19*df0 - 5*df1 + df2)
        y_update = numpy.reshape(y_update, (state_size,1))
        y_state = numpy.append(y_state, y_update, axis=1)

        x_update = x_state
        y_update = y_state

    return [x_state, y_state]
