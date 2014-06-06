'''Bowdtich/Jorgensen discreteness test for punctured torus groups'''
#cython: cdivision=True

from gsl_complex cimport *
from libc.math cimport fabs

#cdef extern from "math.h":
#    double fabs(double x) nogil

cdef int MAX_STEP_TO_SINK=10000
cdef int MAX_DEPTH=30000
cdef double BOWDITCH_THRESHOLD=0.0001
cdef double JORGENSEN_THRESHOLD=0.9999
cdef double HUGE_TRACE=1.0e11

# Status codes for best guess about representation's properties
# Discrete really means "In Bowditch set, did not see any Jorgensen violations along the way"
# Indiscrete means "Found a violation of Jorgensen's inequality, therefore not discrete"
cdef enum rep_type:
    cREP_DISCRETE = 1
    cREP_INDISCRETE = 0
    cREP_UNCERTAIN = 128

REP_DISCRETE = 1
REP_INDISCRETE = 0
REP_UNCERTAIN = 128

cdef double bowditch_H(gsl_complex x) nogil:
    cdef double labs

    labs = gsl_complex_abs(
        gsl_complex_mul_real(
            gsl_complex_add(x,
                            gsl_complex_sqrt (gsl_complex_add_real
                                              (gsl_complex_mul (x, x), -4.0))),
            0.5))
    return fabs(2.0 * (labs + 1.0) / ((labs - 1.0)))


cdef int edge_in_tree(gsl_complex x, gsl_complex y) nogil:
    if (gsl_complex_abs (x) <= (BOWDITCH_THRESHOLD + 3.0)) and \
       (gsl_complex_abs (y) <= (BOWDITCH_THRESHOLD + 1.0 + bowditch_H (x))):
        return 1

    if (gsl_complex_abs (y) <= (BOWDITCH_THRESHOLD + 3.0)) and \
        (gsl_complex_abs (x) <= (BOWDITCH_THRESHOLD + 1.0 + bowditch_H (y))):
         return 1

    return 0

#----------------------------------------------------------------------#
# mtdiscreterec(x,y,z,nodecount,depth)                                 #
#                                                                      #
# This is the recursive core of the discreteness algorithm.  It looks  #
# through the oriented Farey graph corresponding to (x,y,z).  If the   #
# Jorgensen inequality is violated at any stage, the group is judged   #
# indiscrete.                                                          #
#                                                                      #
# RETURN: status code indicating discrete, indiscrete, or failure      #
#         if not NULL, number of nodes visited returned in *nodecount  #
#----------------------------------------------------------------------#

cdef rep_type mtdiscreterec (gsl_complex x, gsl_complex y, gsl_complex z, int *nodecount, int depth) nogil:
    cdef int edge_added
    cdef rep_type r1, r2

    if depth > MAX_DEPTH: 
       # Fail due to max depth
        return cREP_UNCERTAIN

    if gsl_complex_abs (x) < JORGENSEN_THRESHOLD:
        # Jorgensen violation => indiscrete
        return cREP_INDISCRETE

    if (gsl_complex_abs(x) + gsl_complex_abs(y) + gsl_complex_abs(z)) > HUGE_TRACE:
        return cREP_UNCERTAIN

    if nodecount != NULL:
        nodecount[0] = nodecount[0] + 1

    # For each direction not yet explored, call mtdiscreterec to
    # investigate that subtree.  If definite indiscreteness is 
    # indicated, return immediately.  Otherwise keep searching.

    r1 = cREP_DISCRETE
    r2 = cREP_DISCRETE

    if edge_in_tree(x,z):
        r1 = mtdiscreterec (gsl_complex_sub(gsl_complex_mul(x, z), y), z, x,nodecount,depth+1)
        if r1 == cREP_INDISCRETE:
            return r1

    if edge_in_tree(x,y):
        r2 = mtdiscreterec(gsl_complex_sub(gsl_complex_mul(x, y), z), x, y,nodecount,depth+1)
        if r2 == cREP_INDISCRETE:
            return r2

    if (r1 == cREP_UNCERTAIN) or (r2 == cREP_UNCERTAIN):
        # Were undecided before, and no certificate of indiscreteness found later.
        return cREP_UNCERTAIN

    return cREP_DISCRETE

#----------------------------------------------------------------------#
# mtdiscrete(x,y,z,nodecount)                                          #
#                                                                      #
# External interface to the discreteness algorithm.  First look for a  #
# sink, and begin investigating the Farey tree from there using the    #
# recursive procedure mtdiscreterec.                                   #
#                                                                      #
# RETURN: status code indicating discrete, indiscrete, or failure      #
#         if not NULL, number of nodes visited returned in *nodecount  #
#----------------------------------------------------------------------#

cdef rep_type mtdiscrete(gsl_complex x, gsl_complex y, gsl_complex z, int *nodecount) nogil:
    cdef gsl_complex w, t
    cdef rep_type r1, r2, r3
    cdef int n = 0
    cdef int sink = 0

    if nodecount != NULL:
        nodecount[0] = 0

    # ----------------------------------------------------------------------
    # 1. LOCATE SINK
    # ----------------------------------------------------------------------

    w = gsl_complex_rect (3.0, 0.0)
    while (n < MAX_STEP_TO_SINK) and (gsl_complex_abs(w) > JORGENSEN_THRESHOLD) and (sink < 3):
        n = n + 1
        w = gsl_complex_sub( gsl_complex_mul(y, z), x)
        if gsl_complex_abs(w) < gsl_complex_abs(x):
            x = w
            sink = 0
        else:
            sink = sink + 1
        t = x
        x = y
        y = z
        z = t

    if nodecount != NULL:
        nodecount[0] = nodecount[0] + n

    if gsl_complex_abs(w) < JORGENSEN_THRESHOLD:
        return cREP_INDISCRETE

    if n >= MAX_STEP_TO_SINK:
        # Failure to find sink is fatal.  Give up.
        return cREP_UNCERTAIN

    # At this point we know a sink was found successfully.

    # ----------------------------------------------------------------------
    # 2. LOCAL TESTS
    # ----------------------------------------------------------------------

    if (gsl_complex_abs (x) < JORGENSEN_THRESHOLD) or \
       (gsl_complex_abs (y) < JORGENSEN_THRESHOLD) or \
       (gsl_complex_abs (z) < JORGENSEN_THRESHOLD):
        return cREP_INDISCRETE


    # ----------------------------------------------------------------------
    # 3. RECURSIVE CALLS
    # ----------------------------------------------------------------------

    # If no edges are examined, we will fall back to this default return value
    r1 = cREP_DISCRETE
    r2 = cREP_DISCRETE
    r3 = cREP_DISCRETE

    if edge_in_tree(x,y):
        r1 = mtdiscreterec(gsl_complex_sub(gsl_complex_mul(x, y), z), x, y, nodecount, 1)
        if r1 == cREP_INDISCRETE:
            return r1

    if edge_in_tree(x,z):
        r2 = mtdiscreterec(gsl_complex_sub(gsl_complex_mul(x, z), y), z, x, nodecount, 1)
        if r2 == cREP_INDISCRETE:
            return r2

    if edge_in_tree(y,z):
        r3 = mtdiscreterec (gsl_complex_sub(gsl_complex_mul(y, z), x), y, z, nodecount, 1)
        if r3 == cREP_INDISCRETE:
            return r3

    if (r1 == cREP_UNCERTAIN) or (r2 == cREP_UNCERTAIN) or (r3 == cREP_UNCERTAIN):
        return cREP_UNCERTAIN
    
    return cREP_DISCRETE


def classify(*args):
    '''Classify a representation as discrete, indiscrete, or undecided
    using Bowditch+Jorgensen test.

    x,y,z --- Markov triple representing tr(a), tr(b), tr(ab) for
    generators (a,b) of a punctured torus group.

    The arguments must correspond to a type-preserving representation
    where tr([a,b]) = -2.  Equivalently x,y,z must satisfy:

    x**2 + y**2 + z**2 == x*y*z

    Returns REP_DISCRETE, REP_INDISCRETE, or REP_UNCERTAIN
    '''

    try:
        x,y,z = [ complex(t) for t in args ]
    except TypeError:
        x,y,z = [ complex(t) for t in args[0] ]
        
    cdef int nodecount
    cdef rep_type res
    cdef gsl_complex cx, cy, cz

    cx.dat[0] = x.real
    cx.dat[1] = x.imag

    cy.dat[0] = y.real
    cy.dat[1] = y.imag

    cz.dat[0] = z.real
    cz.dat[1] = z.imag

    res = mtdiscrete(cx,cy,cz,&nodecount)
    return res, nodecount

