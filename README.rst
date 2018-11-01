******************************************
CP1: Complex projective structures toolkit
******************************************

Introduction
============

:Author: David Dumas (david@dumas.io, http://dumas.io/)
:Web site: http://github.com/daviddumas/cp1
:Version: 0.0.2
:License: CP1 is released under the GNU General Public License v3 (http://www.gnu.org/licenses/gpl-3.0.txt)

Purpose
=======

CP1 is a Python module for computing holonomy representations of
complex projective structures on Riemann surfaces, and for testing
these representations for discreteness.

CP1 is intended to replicate some of the core functionality of the
software package Bear (http://bear.sf.net/), making it easier to use
these functions in Python programs.


Development status
==================

This is beta code.


Installation
============

Dependencies
------------

* Python 3.6+ or 2.6+
* Numpy (http://www.numpy.org/)
* GCC
* GSL, the GNU Scientific Library (http://www.gnu.org/software/gsl)
* Cython (optional, only necessary if you want to edit the ``.pyx`` and ``.pxd``
  source files)

Installing CP1
--------------

On a recent Debian or Ubuntu linux distribution, these commands will
ensure the dependencies are available and then install the module
system-wide (run them in the root directory of the CP1 source code):

::

    sudo apt-get install build-essential python-dev python-numpy libgsl0-dev
    sudo python setup.py install


Name
====

The package name is CP1, which represents the mathematical notation
for the complex projective line.  (The number 1 should be a
superscript.)  The "tarname", or name of the package as rendered in
filenames and for importing the python module, is ``cp1``.


Acknowledgements
================

CP1 uses Cython wrappers for the ``gsl_odeiv.h`` and ``gsl_complex.h``
GSL headers that were adapted from the CythonGSL project by Thomas Wiecki
(http://github.com/twiecki/CythonGSL), which is itself a fork of the PyrexGsl by Mario Pernici
(http://wwwteor.mi.infn.it/~pernici/pyrexgsl/pyrexgsl.html).
