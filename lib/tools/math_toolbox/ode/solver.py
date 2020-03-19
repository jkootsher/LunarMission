import lib.tools.math_toolbox.integration.integrands as seed
import lib.tools.math_toolbox.integration.integrals as integrate

def ode113v(**kwargs):
    ''' Adam's & Bashforth ODE (non-stiff) for single state dynamics '''
    kwargs['h'] = kwargs['step_size']
    kwargs['x'] = kwargs['samples']
    kwargs['y'] = kwargs['state']
    return integrate.ABM4(seed.CRTBP.vfhandle, **kwargs)

def ode113m(**kwargs):
    ''' Adam's & Bashforth ODE (non-stiff) for multi-state dynamics '''
    kwargs['h'] = kwargs['step_size']
    kwargs['x'] = kwargs['samples']
    kwargs['y'] = kwargs['state']
    return integrate.ABM4(seed.CRTBP.mfhandle, **kwargs)

def test(**kwargs):
    ''' Verify numerical integration '''
    kwargs['h'] = kwargs['step_size']
    kwargs['x'] = kwargs['samples']
    kwargs['y'] = kwargs['state']
    return integrate.RK4(seed.CRTBP.testhandle, **kwargs)