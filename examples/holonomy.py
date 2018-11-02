#!/usr/bin/env python
'''Compute holonomy of a complex projective structure'''

from __future__ import print_function
import os
import sys
import argparse

# Support module built in-place (in ../cp1) or system-wide
try:
    import cp1
except ImportError:
    sys.path.insert(0, '..')
    import cp1

parser = argparse.ArgumentParser(description='Compute holonomy of a projective structure')
group = parser.add_mutually_exclusive_group()
group.add_argument('--s04',action='store_true',
                   help='Surface = four puncture sphere')
group.add_argument('--t11',action='store_true',
                   help='Surface = punctured torus')

group2 = parser.add_mutually_exclusive_group()
group2.add_argument('-L','--lambda',type=complex,dest='L',
                   help='Point in moduli space C - {0,1}')
group2.add_argument('-t','--tau',type=complex,
                   help='Point in Teichmuller space (upper half plane)')

parser.add_argument('-c','--coef',type=complex,dest='C',
                    help='Complex coefficient determining the projective structure')

args = parser.parse_args()

if not args.s04 and not args.t11:
    print('Warning: Topological type not specified. Using four-punctured sphere (S04).\n')
    args.s04 = True

if args.tau == None and args.L == None:
    print('Warning: Riemann surface not specified.  Using tau = lambda = exp(pi*i/3).\n')
    args.L = 0.5+0.866025j

if args.C == None:
    print('Warning: Coefficient not specified.  Using C = -i/sqrt(3).\n')
    args.C = -0.57735j


if args.s04:
    if args.L:
        m = cp1.s04_lambda_hol(L=args.L,C=args.C,contours=3)
    else:
        m = cp1.s04_tau_hol(tau=args.tau,C=args.C,contours=3)
    m_t11 = cp1.s04_to_t11(m)
    d,n = cp1.classify(m_t11)
else:
    if args.L:
        m = cp1.t11_lambda_hol(L=args.L,C=args.C,contours=3)
    else:
        m = cp1.t11_tau_hol(tau=args.tau,C=args.C,contours=3)
    d,n = cp1.classify(m)


if args.s04:
    print('Projective structure on S04:')
else:
    print('Projective structure on T11:')

if args.tau:
    print('  tau =\t',args.tau)
    print('  lam =\t',cp1.modularlambda(args.tau))
else:
    print('  lam =\t',args.L)

print('  C =\t',args.C)
print('Holonomy traces:')
print('  x =\t',m[0])
print('  y =\t',m[1])
print('  z =\t',m[2])
print('Holonomy type:')

if d == cp1.REP_DISCRETE:
    print('  Probably DISCRETE (satisfies Bowdtich condition, no Jorgensen violation found)')
    print('  Examined %d Farey triangles' % n)
elif d == cp1.REP_INDISCRETE:
    print('  NOT DISCRETE (violation of Jorgensen inequality was found)')
    print('  Examined %d Farey triangles' % n)
else:
    print('  Unable to determine')
    print('  Examined %d Farey triangles' % n)


