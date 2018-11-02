#!/usr/bin/env python

'''Create bitmap image of an extended Bers slice showing discrete and non-discrete holonomy.'''

from __future__ import print_function
import os
import sys

def cplx_pp(z):
    s = '%g' % z.real
    if z.imag > 0:
        s += ' + %g i' % z.imag
    else:
        s += ' - %g i' % (-z.imag)
    return s

# Support module built in-place (in ../cp1) or system-wide
try:
    from cp1 import (modularlambda, t11_tau_hol, classify, HolonomyException,
                     REP_DISCRETE, REP_INDISCRETE, REP_UNCERTAIN)
except ImportError:
    sys.path.insert(0, '..')
    from cp1 import (modularlambda, t11_tau_hol, classify, HolonomyException,
                     REP_DISCRETE, REP_INDISCRETE, REP_UNCERTAIN)

try:
    from PIL import Image, ImageDraw
except ImportError:
    sys.stderr.write('''Unable to load PIL.

This example requires the Python Imaging Library (PIL).
See:  http://www.pythonware.com/products/pil/

''')
    sys.exit(1)

import argparse

parser = argparse.ArgumentParser(description='Draw a Bers slice image for a point in Teichmuller space')
parser.add_argument('-t','--tau',type=complex,default=0.5+0.866j,
                    help='Point in Teichmuller space (upper half plane)')
parser.add_argument('-c','--center',type=complex,default=-0.57735026918962576451j,
                    help='Center of the image (complex number giving a quadratic differential)')
parser.add_argument('-r','--radius',type=float,default=0.5,
                    help='Radius of ths image')
parser.add_argument('-s','--size',type=int,default=256,
                    help='Size of the image (NxN pixels, specify N)')
parser.add_argument('-o','--output',
                    help='Output filename (default is based on current time)')

args = parser.parse_args()

if args.output is None:
    from time import time
    args.output = 'bs%d.png' % int(time())
    print('Using output filename "%s"' % args.output)

img = Image.new('RGB',(args.size, args.size + 70), (255,255,255))

L = modularlambda(args.tau)

# Say what we are doing
print('Creating %dx%d Bers slice image:' % (args.size,args.size))
print('\ttau = %s' % cplx_pp(args.tau))
print('\tlambda = %s' % cplx_pp(L))
print('\tcenter = %s' % cplx_pp(args.center))
print('\tradius = %g' % args.radius)

# Prepare the label text
lines = []
lines.append('tau = %s' % cplx_pp(args.tau))
lines.append('lambda = %s' % cplx_pp(modularlambda(args.tau)))
lines.append('ctr = %s' % cplx_pp(args.center))
lines.append('radius = %g' % args.radius)
d = ImageDraw.Draw(img)
xpos = 5
ypos = args.size + 5
for line in lines:
    w,h = d.textsize(line)
    d.text((xpos,ypos),line,fill=(0,0,0))
    ypos = ypos + h + 3

d.rectangle( [(0,0),(args.size-1,args.size-1)], fill=(180,180,180) )

# Draw the slice
pix = img.load()
try:
    for j in range(args.size):
        print('Row %d of %d' % (j+1,args.size))
        for i in range(args.size):
            tx = 2.0*(float(i)/(args.size-1)) - 1.0
            ty = 1.0 - 2.0*(float(j)/(args.size-1))
            C = args.center + tx*args.radius + 1j*ty*args.radius
            try:
                mt = t11_tau_hol(tau=args.tau, C=C)
                d, n = classify(mt)
                if d is REP_DISCRETE:
                    pix[i,j] = (0,0,0)
                elif d is REP_INDISCRETE:
                    k = n
                    if k > 200:
                        k = 200
                    pix[i,j] = (255-k,255-k,255)
                else:
                    pix[i,j] = (255,0,0)
            except HolonomyException:
                pix[i,j] = (0,255,0)
except KeyboardInterrupt:
    print('Interrupted: Will try to write partial progress to file.')
    pass

# Write to a file
print('Writing "%s"' % args.output)
img.save(args.output)

        
        
