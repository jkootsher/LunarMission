from lib.general.datatypes import Vector

def jacobi(r=Vector(), v=Vector(), mu=1):
    ''' Calculate the Jacobi energy for the current orbit state '''
    mu_vector = Vector((mu,0.0,0.0), r.frame)
    r1_vector = r + mu_vector
    r1_norm = r1_vector.norm()

    mu_vector = Vector((mu-1,0.0,0.0), r.frame)
    r2_vector = r + mu_vector
    r2_norm = r2_vector.norm()

    u = -.5*v.norm()**2 + 2*(r.i.ts[r.n]**2 + r.j.ts[r.n]**2 + (1-mu)/r1_norm + mu/r2_norm)
    return u