'''Integrate a projective connection on a four-punctured sphere'''

from gsl_complex cimport *
from gsl_odeiv cimport *

ctypedef struct pc_seg_odef_params:
    gsl_complex L
    gsl_complex C
    gsl_complex p0
    gsl_complex p1

# Necessary to include "const" qualifier in generated c
cdef extern from *:
    ctypedef double const_double "const double"

# Const in this prototype eliminates a compiler warning
cdef int pc_seg_odef (double t, const_double y[], double f[], void *P) nogil:
    # cast parameters
    cdef pc_seg_odef_params *params
    params = (<pc_seg_odef_params *>P)

    cdef gsl_complex L,C,p0,p1
    L = params.L
    C = params.C
    p0 = params.p0
    p1 = params.p1
    
    # unpack y-values into complex u-values
    cdef gsl_complex u[4]
    u[0] = gsl_complex_rect(y[0],y[1])
    u[1] = gsl_complex_rect(y[2],y[3])
    u[2] = gsl_complex_rect(y[4],y[5])
    u[3] = gsl_complex_rect(y[6],y[7])

    # compute matrix coefficients
    cdef gsl_complex a, z

    a = gsl_complex_sub(p1,p0)
    z = gsl_complex_add(p0,gsl_complex_mul_real(a,t))

    cdef gsl_complex z1,zL
    z1 = gsl_complex_sub_real(z,1.0)
    zL = gsl_complex_sub(z,L);

    cdef gsl_complex q1,q2,q3
    q1 = gsl_complex_inverse(gsl_complex_mul(gsl_complex_mul(z,z),gsl_complex_mul(z1,z1)))
    q2 = gsl_complex_inverse(gsl_complex_mul(zL,zL))
    q3 = gsl_complex_div(C,gsl_complex_mul(z,gsl_complex_mul(z1,zL)))
    
    cdef gsl_complex phi
    phi = gsl_complex_mul_real(gsl_complex_add(gsl_complex_mul_real(gsl_complex_add(q1,q2),0.5),q3),-0.5)

    # phi = -0.25/((z*(z-1))**2) - 0.25/(z-L)**2 - 0.5*C/(z*(z-1)*(z-L))

    # store in complex uprime-values
    cdef gsl_complex uprime[4]
    uprime[0] = gsl_complex_mul(a,u[1])
    uprime[1] = gsl_complex_mul(a,gsl_complex_mul(phi,u[0]))
    uprime[2] = gsl_complex_mul(a,u[3])
    uprime[3] = gsl_complex_mul(a,gsl_complex_mul(phi,u[2]))

    # pack into double return values
    f[0] = GSL_REAL(uprime[0])
    f[1] = GSL_IMAG(uprime[0])
    f[2] = GSL_REAL(uprime[1])
    f[3] = GSL_IMAG(uprime[1])
    f[4] = GSL_REAL(uprime[2])
    f[5] = GSL_IMAG(uprime[2])
    f[6] = GSL_REAL(uprime[3])
    f[7] = GSL_IMAG(uprime[3])

    return GSL_SUCCESS


cdef void make_eye(double *m) nogil:
    m[0] = 1.0  # re(a)
    m[1] = 0.0  # im(a)
    m[2] = 0.0  # re(b)
    m[3] = 0.0  # im(b)
    m[4] = 0.0  # re(c)
    m[5] = 0.0  # im(c)
    m[6] = 1.0  # re(d)
    m[7] = 0.0  # im(d)

cdef class pcint:
    '''Integrate a projective connection along a segment, for the 4-punctured sphere universal family'''
    cdef gsl_odeiv_step *step
    cdef gsl_odeiv_control *control
    cdef gsl_odeiv_evolve *evolve
    cdef gsl_odeiv_system sys
    cdef pc_seg_odef_params P
    cdef double initstep
    cdef int maxstep

    def __cinit__(self, double atol=0.00001, double rtol=0.0, double initstep=0.01, int maxstep=5000):
        self.initstep = initstep
        self.maxstep = maxstep
        self.step = gsl_odeiv_step_alloc(gsl_odeiv_step_rk8pd, 8)  # (type, rank)
        self.control = gsl_odeiv_control_standard_new(atol, rtol, 0.0, 0.0)
        self.evolve  = gsl_odeiv_evolve_alloc(8)
        self.sys.function = pc_seg_odef
        self.sys.jacobian = NULL
        self.sys.dimension = 8
        self.sys.params = &self.P
    
    def __dealloc__(self):
        if self.step != NULL:
            gsl_odeiv_step_free(self.step)
        if self.control != NULL:
            gsl_odeiv_control_free(self.control)
        if self.evolve != NULL:
            gsl_odeiv_evolve_free(self.evolve)

    def seg_int(self, L, C, p0, p1, init=None, maxstep=None):
        self.P.L = gsl_complex_rect(L.real, L.imag)
        self.P.C = gsl_complex_rect(C.real, C.imag)
        self.P.p0 = gsl_complex_rect(p0.real, p0.imag)
        self.P.p1 = gsl_complex_rect(p1.real, p1.imag)

        cdef double t, h, y[8]
        t = 0.0
        h = self.initstep

        if init == None:
            make_eye(y)
        else:
            y[0] = init[0][0].real
            y[1] = init[0][0].imag
            y[2] = init[0][1].real
            y[3] = init[0][1].imag
            y[4] = init[1][0].real
            y[5] = init[1][0].imag
            y[6] = init[1][1].real
            y[7] = init[1][1].imag

        cdef int MS
        if maxstep == None:
            MS = self.maxstep
        else:
            MS = maxstep

        cdef int status
        cdef int n
        n = 0
        status = GSL_SUCCESS
        while (t < 1.0) and (n < MS):
            status = gsl_odeiv_evolve_apply(self.evolve,
                                            self.control,
                                            self.step,
                                            &self.sys,
                                            &t, 1.0, &h, y)

            if (status != GSL_SUCCESS):
                break

            n = n + 1

        if (n == MS) or (status != GSL_SUCCESS):
            return None

        return [ [ y[0] + 1j*y[1], y[2] + 1j*y[3] ],
                 [ y[4] + 1j*y[5], y[6] + 1j*y[7] ] ]


LAMBDA_ITER_MAX = 200
DBL_EPSILON = 2.2204460492503131e-16
M_PI = 3.14159265358979323846264338328
M_PI_4 = 0.78539816339744830961566084582

cdef gsl_complex _modularlambda(gsl_complex tau):
    cdef int n=0

    cdef gsl_complex q = gsl_complex_exp(gsl_complex_mul_imag(tau,M_PI))
    cdef gsl_complex q2 = gsl_complex_mul(q,q)
    cdef gsl_complex q14 = gsl_complex_exp(gsl_complex_mul_imag(tau,M_PI_4))

    cdef gsl_complex accum = gsl_complex_rect(0.0,0.0)
    cdef gsl_complex nextm = q2
    cdef gsl_complex qpower = gsl_complex_rect(1.0,0.0)

    while ((gsl_complex_abs(qpower) > 2.0*DBL_EPSILON) and (n < LAMBDA_ITER_MAX)):
        accum = gsl_complex_add(accum, qpower)
        qpower = gsl_complex_mul(qpower, nextm)
        nextm = gsl_complex_mul(nextm, q2)
        n = n + 1

    if (n >= LAMBDA_ITER_MAX):
        return gsl_complex_rect(0.0,0.0)

    cdef gsl_complex theta2 = gsl_complex_mul_real(gsl_complex_mul(q14,accum),2.0)
  
    accum = gsl_complex_rect(0.5,0.0)
    nextm = gsl_complex_mul(q,q2)
    qpower = q
  
    n = 0
    while ((gsl_complex_abs(qpower) > 2.0*DBL_EPSILON) and (n < LAMBDA_ITER_MAX)):
        accum = gsl_complex_add(accum, qpower)
        qpower = gsl_complex_mul(qpower, nextm)
        nextm = gsl_complex_mul(nextm, q2)
        n = n + 1

    if (n >= LAMBDA_ITER_MAX):
        return(gsl_complex_rect(0.0,0.0));

    cdef gsl_complex theta3 = gsl_complex_mul_real(accum,2.0)

    cdef gsl_complex x = gsl_complex_div(theta2,theta3)
    return gsl_complex_mul(gsl_complex_mul(x,x), gsl_complex_mul(x,x))

def modularlambda(t):
    cdef gsl_complex tau, L
    tau = gsl_complex_rect(t.real,t.imag)
    L = _modularlambda(tau)
    return complex(L.dat[0],L.dat[1])
