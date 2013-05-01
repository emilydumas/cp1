import cp1

TESTDELTA = 0.00001

def test_t11_hex():
    L = 0.5 + 0.86602540378443864676j
    C = -0.57735026918962576451j
    xref = 3.0
    yref = 3.0
    zref = 3.0
    x,y,z = cp1.t11_lambda_hol(L,C,contours=3)
    assert( abs(x-xref) < TESTDELTA)
    assert( abs(y-yref) < TESTDELTA)
    assert( abs(z-zref) < TESTDELTA)

def test_t11_square():
    L = 0.5 
    C = 0.0
    xref = 2.82842712474619009760
    yref = 2.82842712474619009760
    zref = 4.0
    x,y,z = cp1.t11_lambda_hol(L,C,contours=3)
    assert( abs(x-xref) < TESTDELTA)
    assert( abs(y-yref) < TESTDELTA)
    assert( abs(z-zref) < TESTDELTA)
