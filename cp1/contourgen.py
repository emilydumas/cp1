'''Integration contours for the four-puncture sphere'''

# ----------------------------------------------------------------------
# SIMPLE CONTOURS
#
# 4,4,6 vertices, work for lambda in a rectangle 
# [epsilon,1-epsilon] x [-big,big]
#
# Will fail when lambda is very close to imaginary axis
# ----------------------------------------------------------------------

def _make_rect(xmin,xmax,ymin,ymax):
    return [ xmax+1j*ymax,
             xmin+1j*ymax,
             xmin+1j*ymin,
             xmax+1j*ymin ]

def simp1(L):
    # Rectangular loop around 0 and L ( = lambda )
    xmin = -0.25
    xmax = 0.5*(L.real + 1.0)
    if (L.real + 0.25) < xmax:
        xmax = L.real + 0.25
    ymin = -0.25;
    if (L.imag - 0.25) < ymin:
        ymin = L.imag - 0.25
    ymax = 0.25;
    if (L.imag + 0.25) > ymax:
        ymax = L.imag + 0.25
    return _make_rect(xmin,xmax,ymin,ymax)

def simp2(L):
    # Rectangular loop around L and 1
    xmin = 0.5*L.real
    xmax = 1.25
    if (L.real - 0.25) > xmin:
        xmin = L.real - 0.25
    ymin = -0.25
    if (L.imag - 0.25) < ymin:
        ymin = L.imag - 0.25
    ymax = 0.25
    if (L.imag + 0.25) > ymax:
        ymax = L.imag + 0.25
    return _make_rect(xmin,xmax,ymin,ymax)

def simp3(L):
    # Hexagonal loop with 0, 1 inside, L outside
    xmin = -0.25
    xmax = 1.25
    xmid1 = 0.5*L.real
    if (L.real - 0.25) > xmid1:
        xmid1 = L.real - 0.25
    xmid2 = 0.5*(L.real + 1.0)
    if (L.real + 0.25) < xmid2:
        xmid2 = L.real + 0.25
    ymin = -0.25
    if (L.imag - 0.25) < ymin:
        ymin = L.imag - 0.25
    ymax = 0.25
    if (L.imag + 0.25) > ymax:
        ymax = L.imag + 0.25
    return [xmax+1j*ymax,
            xmid2+1j*ymax,
            xmid2+1j*ymin,
            xmid1+1j*ymin,
            xmid1+1j*ymax,
            xmin+1j*ymax,
            xmin+1j*ymin,
            xmax+1j*ymin]

def simple(L,n):
    return [ gamma(L) for gamma in [simp1,simp2,simp3][:n] ]

# ----------------------------------------------------------------------
# ADVANCED CONTOURS
#
# more vertices, more logic, but work for lambda in 
# 
# ( [-0.2,0.75] x [-big,big] ) - [-0.2,0]
#
# Note that this includes a fundamental domain for the modular group,
# but excludes some values that work with the simple contours above
# ----------------------------------------------------------------------

def _signum(x):
    if x<0.0:
        return -1
    elif x> 0.0:
        return 1
    return 0

def horiz_int(L,y):
    if abs(L) > 100.0*abs(L.real):
        return -100.0*_signum(L.imag * y)
    return (-2.0*y*L.imag + abs(L)**2)/(2.0*L.real)

def vert_int(L,x):
    return (abs(L)**2 - 2.0*x*L.real)/(2.0*L.imag)

def advanced(L,n):
    #print 'L = ',L
    contours = [simp1(L)]

    M = L
    if L.real < 1e-8:
        M = L.imag * 1j

    ytop = max(0.5,L.imag+0.25)
    ybot = min(-0.5,L.imag-0.25)
    ixtop = horiz_int(M,ytop)
    ixbot = horiz_int(M,ybot)

    #print ytop, ybot, ixtop, ixbot

    gamma = [1.25 + ytop*1j, 1.25 + ybot*1j]
    if ixbot < -0.25:
        gamma.extend([-0.25 + ybot*1j,-0.25 + vert_int(M,-0.25)*1j])
    elif ixbot > 0.5:
        gamma.extend([0.5 + ybot*1j, 0.5 + vert_int(M,0.5)*1j])
    else:
        gamma.extend([ixbot + ybot*1j])

    if ixtop < -0.25:
        gamma.extend([-0.25 + vert_int(M,-0.25)*1j,-0.25 + ytop*1j])
    elif ixtop > 0.5:
        gamma.extend([0.5 + vert_int(M,0.5)*1j,0.5 + ytop*1j])
    else:
        gamma.extend([ixtop + ytop*1j])

    #print 'gamma_2 = ',gamma
    contours.append(gamma)

    if n < 3:
        return contours

    ybot = min(-0.25,L.imag-0.25)
    ixbot = horiz_int(M,ybot)
    gamma = [-0.25 + ybot*1j, 1.25 + ybot*1j, 1.25 + ytop*1j, 0.75 + ytop*1j]

    if ixbot < -0.2:
        gamma.extend([0.75 + ybot*1j,-0.2 + ybot*1j,-0.2 + vert_int(M,-0.2)*1j])
    elif ixbot > 0.75:
        gamma.extend([0.75 + vert_int(M,0.75)*1j])
    else:
        gamma.extend([0.75 + ybot*1j, ixbot + ybot*1j])

    if ixtop < -0.25:
        gamma.extend([-0.25 + vert_int(M,-0.25)*1j])
    elif ixtop > 0.5:
        gamma.extend([0.5 + vert_int(M,0.5)*1j, 0.5 + ytop*1j,-0.25 + ytop*1j])
    else:
        gamma.extend([ixtop + ytop*1j, -0.25 + ytop*1j])

    #print 'gamma_2 = ',gamma
    contours.append(gamma)

    return contours
