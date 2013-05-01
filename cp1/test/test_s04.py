import cp1
import pytest

TESTDELTA = 0.00001

@pytest.mark.parametrize( ("fund_domain_contours"), [True, False])
def test_s04_hex(fund_domain_contours):
    L = 0.5 + 0.86602540378443864676j
    C = -0.57735026918962576451j
    xref = -7.0
    yref = -7.0
    zref = -7.0
    x,y,z = cp1.s04_lambda_hol(L,C,contours=3,tol=1e-8,
                               fund_domain_contours=fund_domain_contours)
    assert( abs(x-xref) < TESTDELTA)
    assert( abs(y-yref) < TESTDELTA)
    assert( abs(z-zref) < TESTDELTA)

@pytest.mark.parametrize( ("fund_domain_contours"), [True, False])
def test_s04_square(fund_domain_contours):
    L = 0.5
    C = 0.0
    xref = -6.0
    yref = -6.0
    zref = -14.0
    x,y,z = cp1.s04_lambda_hol(L,C,contours=3,tol=1e-8,
                               fund_domain_contours=fund_domain_contours)
    assert( abs(x-xref) < TESTDELTA)
    assert( abs(y-yref) < TESTDELTA)
    assert( abs(z-zref) < TESTDELTA)

@pytest.mark.parametrize( ("fund_domain_contours"), [True, False])
def test_s04_largeC(fund_domain_contours):
    L = 0.5 + 0.86602540378443864676j
    C = 2.0 + 3.0j - L
    xref = -222.533070-262.929905j
    yref = 1272.251270-2334.830978j
    zref = -2.848143-7.056932j
    x,y,z = cp1.s04_lambda_hol(L,C,contours=3,tol=1e-8,
                               fund_domain_contours=fund_domain_contours)
    assert( abs(x-xref) < TESTDELTA)
    assert( abs(y-yref) < TESTDELTA)
    assert( abs(z-zref) < TESTDELTA)

@pytest.mark.parametrize( ("fund_domain_contours"), [True, False])
def test_s04_extreme(fund_domain_contours):
    L = 0.00001 + 0.00001j
    C = 0.5 + 0.3j - L
    xref = 0.861380+5.348997j
    yref = -8013.211655-1721.651316j
    zref = 519.186599-1373.061159j
    x,y,z = cp1.s04_lambda_hol(L,C,contours=3,tol=1e-8,
                               fund_domain_contours=fund_domain_contours)
    assert( abs(x-xref) < TESTDELTA)
    assert( abs(y-yref) < TESTDELTA)
    assert( abs(z-zref) < TESTDELTA)

def test_s04_compare_simp_adv():
    import cmath
    N = 20
    for i in range(N):
        t = float(i)/float(N-1)
        L = 0.3*(1.0-cmath.exp(2.0j*cmath.pi*(0.01 + 0.98*t)))
        C = 1.2+0.5j+cmath.exp(4.0j*cmath.pi*t)
        x1,y1,z1 = cp1.s04_lambda_hol(L,C,contours=3,tol=1e-8,
                                      fund_domain_contours=False)
        x2,y2,z2 = cp1.s04_lambda_hol(L,C,contours=3,tol=1e-8,
                                      fund_domain_contours=False)
        assert( abs(x1-x2) < TESTDELTA)
        assert( abs(y1-y2) < TESTDELTA)
        assert( abs(z1-z2) < TESTDELTA)


    
