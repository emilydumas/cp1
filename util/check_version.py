#!/usr/bin/python

'''Check that version numbers are consistent'''

import os
import re

cp1root = os.path.realpath(os.path.dirname(os.path.realpath(__file__))+'/..')

# setup.py
setup_version_re = re.compile('\s*version\s*=\s*\'([^\']+)\'')
f = file(cp1root + '/setup.py','rt')
for line in f:
    m = setup_version_re.search(line)
    if m:
        print 'setup.py:\t%s' % m.group(1)
        break  
else:
    print 'ERROR: No version found in setup.py'
  
# README.rst
readme_re = re.compile(':Version:\s*(\S+)$')
f = file(cp1root + '/README.rst','rt')
for line in f:
    m = readme_re.search(line)
    if m:
        print 'README.rst:\t%s' % m.group(1)
        break
else:
    print 'ERROR: No version found in README.rst'

# __init__.py
init_re = re.compile('__version__\s*=\s*\'([^\']+)\'')
f = file(cp1root + '/cp1/__init__.py','rt')
for line in f:
    m = init_re.search(line)
    if m:
        print '__init__.py:\t%s' % m.group(1)
        break
else:
    print 'ERROR: No version found in __init__.py'

