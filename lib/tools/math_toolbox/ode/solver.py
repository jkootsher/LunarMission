import lib.tools.math_toolbox.integration.integrands as seed
import lib.tools.math_toolbox.integration.integrals as integrate

from lib.simcore.function_library import feval


def ode113v(**kwargs):
    ''' Adam's & Bashforth ODE (non-stiff) for single state dynamics '''
    kwargs['h'] = kwargs['delta']
    kwargs['x'] = kwargs['tspan']
    kwargs['y'] = kwargs['state']
    return integrate.ABM4(seed.CRTBP.vfhandle, **kwargs)

def ode113m(**kwargs):
    ''' Adam's & Bashforth ODE (non-stiff) for multi-state dynamics '''
    kwargs['h'] = kwargs['delta']
    kwargs['x'] = kwargs['tspan']
    kwargs['y'] = kwargs['state']
    return integrate.ABM4(seed.CRTBP.mfhandle, **kwargs)

def newton_raphson(fhandle=None, **kwargs):
    ''' Newton Raphson root finding method '''
    dfhandle = kwargs['df']
    x_initial = kwargs['x']

    tolerance = 1e-8
    MAX_ITERS = 1000

    curr_iter = 0
    f = feval(fhandle, **kwargs)
    df = feval(dfhandle, **kwargs)
    
    x = x_initial - f/df
    error = abs(x - x_initial)

    while error > tolerance and curr_iter < MAX_ITERS:
        kwargs['x'] = x
        f = feval(fhandle, **kwargs)
        df = feval(dfhandle, **kwargs)
        
        x_n = x - f/df
        error = abs(x_n - x)

        curr_iter = curr_iter + 1
        x = x_n
    return x