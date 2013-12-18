'''Holonomy of projective structures on a punctured torus'''
import s04
from numpy import trace
from cmath import sqrt

def markov_z(x,y):
    return 0.5*(x*y - sqrt( (x*y)**2 - 4.0*(x**2 + y**2) ))

def s04_to_t11(traces):
    if len(traces) == 2:
        x,y = [ sqrt(2-t) for t in traces ]
        z = markov_z(x,y)
        return x,y,z
    else:
        # length == 3
        x,y,z = [ sqrt(2-t) for t in traces ]
        if abs(x**2 + y**2 + z**2 + x*y*z) < abs(x**2 + y**2 + z**2 - x*y*z):
            z = -z
        return x,y,z   

# Convention: L = lambda = cross ratio parameter of 4-punctured sphere commensurable
#                          with the desired punctured torus
#             C = quadratic differential parameter
#
# Together pairs (L,C) parameterize all CP^1 structures of bounded
# type on the punctured torus

# This class will compute holonomy for a given L and any number of
# values of C.  (Use this if L will be fixed for many calls, to avoid
# recomputing trajectories each time.)
class T11(s04.S04):
    '''Represents the space of CP1 structures on a punctured torus Riemann surface'''

    def gens(self,C):
        raise NotImplementedError('Matrix generators are not available for T11 holonomy.')

    def traces(self,C):
        '''Compute traces of holonomy group generators'''
        s04traces = [ trace(m) for m in s04.S04.gens(self,C) ]
        return s04_to_t11(s04traces)

# Convenience function to compute trace tuple directly from L and C
# (use this when L and C will change independently between calls)
def t11_lambda_hol(L=0.5,C=0.5,**kwargs):
    '''Compute traces of holonomy group generators for (L,C) projective connection'''
    return s04_to_t11(s04.s04_lambda_hol(L=L,C=C,**kwargs))

def t11_tau_hol(tau=1j,C=0.0,**kwargs):
    '''Compute traces of holonomy group generators for (tau,C) projective connection'''
    return s04_to_t11(s04.s04_tau_hol(tau=tau,C=C,**kwargs))
