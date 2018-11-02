#!/usr/bin/env python
'''Draw ascii art picture of square torus Bers embedding'''

from __future__ import print_function
import os,sys

# Support module built in-place (in ../cp1) or system-wide
try:
    from cp1 import t11_lambda_hol, classify, REP_DISCRETE, REP_INDISCRETE, REP_UNCERTAIN
except ImportError:
    sys.path.insert(0, '..')
    from cp1 import t11_lambda_hol, classify, REP_DISCRETE, REP_INDISCRETE, REP_UNCERTAIN

xsize = 50
ysize = 40

# Square torus is commensurable with C - {0,1,L} for L = 0.5

L = 0.5  

# Accessory parameter is 0 here (center on Fuchsian point)
C0 = 0

# Slightly larger than Bers embedding
radius = 0.16

for j in range(ysize):
    for i in range(xsize):
        tx = 2.0*(float(i)/(xsize-1)) - 1.0
        ty = 1.0 - 2.0*(float(j)/(ysize-1))
        C = C0 + tx*radius + 1j*ty*radius
        mt = t11_lambda_hol(L=0.5, C=C)
        d, n = classify(mt)
        if d == REP_DISCRETE:
            print('*',end='')
        elif d == REP_INDISCRETE:
            print(' ',end='')
        else:
            print('?',end='')
    print('')


        
        
