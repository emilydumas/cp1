from distutils.core import setup,Command
from distutils.extension import Extension

try:
    from Cython.Distutils import build_ext
except ImportError:
    use_cython = False
else:
    use_cython = True

def get_gsl_info():
    import sys,os

    if 'GSL_CONFIG' in os.environ:
        gsl_config_exec = os.getenv('GSL_CONFIG')
    else:
        gsl_config_exec = 'gsl-config'

    gsl_cflags_str = os.popen(gsl_config_exec + ' --cflags').read()
    gsl_libs_str = os.popen(gsl_config_exec + ' --libs').read()

    if (not gsl_cflags_str) or (not gsl_libs_str):
        raise Exception('Unable to find gsl library and include paths (set GSL_CONFIG environment variable to location of gsl-config executable).')

    gsl_inc_dirs = [ f[2:] for f in gsl_cflags_str.split() if f.startswith('-I') ]
    gsl_lib_dirs = [ f[2:] for f in gsl_libs_str.split() if f.startswith('-L') ]
    gsl_libs = [ f[2:] for f in gsl_libs_str.split() if f.startswith('-l') ]
    return {'inc_dirs':gsl_inc_dirs,
            'lib_dirs':gsl_lib_dirs,
            'libs':gsl_libs }

class PyTest(Command):
    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        import sys,subprocess
        errno = subprocess.call([sys.executable, 'runtests.py'])
        raise SystemExit(errno)

cmdclass = { 'test': PyTest }
ext_modules = [ ]

if use_cython:
    extension_sources = { 'pcint': [ 'cp1/pcint.pyx', 'cp1/gsl_odeiv.pxd', 'cp1/gsl_complex.pxd' ],
                          'bowditch': [ 'cp1/bowditch.pyx', 'cp1/gsl_complex.pxd' ] }
    cmdclass.update({ 'build_ext': build_ext })
else:
    extension_sources = { 'pcint': [ 'cp1/pcint.c' ],
                          'bowditch': [ 'cp1/bowditch.c' ] }

gsl_info = get_gsl_info()

ext_modules += [
    Extension('cp1.pcint', extension_sources['pcint'],
              include_dirs=gsl_info['inc_dirs'],
              library_dirs=gsl_info['lib_dirs'],
              libraries=gsl_info['libs']),
    Extension('cp1.bowditch', extension_sources['bowditch'],
              include_dirs=gsl_info['inc_dirs'],
              library_dirs=gsl_info['lib_dirs'],
              libraries=gsl_info['libs']),
    ]

def readme():
    with open('README.rst') as f:
        return f.read()

setup(
    name = 'cp1',
    version = '0.0.1dev',
    description = 'CP1: Complex projective structures toolkit',
    long_description = readme(),
    author = 'David Dumas',
    author_email='david@dumas.io',
    url='http://github.com/daviddumas/cp1',
    license='GPLv3',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Mathematics'
        ],
    cmdclass = cmdclass,
    ext_modules = ext_modules,
    packages = ['cp1']
    )
