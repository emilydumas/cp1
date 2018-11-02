#!/usr/bin/env python
'''Read re(lambda) im(lambda) re(C) im(C) tuples from stdin or file and compute holonomy'''

# This example is mainly for benchmarking holonomy computation

from __future__ import print_function
import os
import sys

# Support module built in-place (in ../cp1) or system-wide
try:
    from cp1 import s04_lambda_hol
except ImportError:
    sys.path.insert(0, '..')
    from cp1 import s04_lambda_hol

# Input filename can be specified; must be last command line argument
if len(sys.argv)>1 and not sys.argv[-1].startswith('-'):
    infile = open(sys.argv[-1],'rt')
else:
    infile = sys.stdin

if '--quiet' in sys.argv:
    quiet = True
else:
    quiet = False

for line in infile:
    fields = [float(x) for x in line[:-1].split()]
    L = fields[0] + 1j*fields[1]
    C = fields[2] + 1j*fields[3]
    t = s04_lambda_hol(L,C)
    if not quiet:
        print('lambda = %f + I %f' % (L.real,L.imag))
        print('C = %f + I %f' % (C.real + L.real,C.imag + L.imag))
        print('pre_mt[0] = %f + I %f' % (t[0].real,t[0].imag))
        print('pre_mt[1] = %f + I %f' % (t[1].real,t[1].imag))
