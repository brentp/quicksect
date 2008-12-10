#from setuptools import setup
#from distutils.extension import Extension
#from Cython.Distutils import build_ext
##
##print 'remember to cython --convert-range intersection.pyx before  setup.py build_ext --inplace'
#import os
#os.system('cython --convert-range src/quicksect.pyx')
#


from distutils.core import setup
from Cython.Distutils.extension import Extension
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(name="quicksect", script_args = ['build',     '--build-base', 'src/',
                     'build_ext', '--force', '--inplace',
                     ],
      ext_modules = [Extension("quicksect", ["src/quicksect.pyx"],
                               include_dirs=[])],
      cmdclass    = {'build_ext' : build_ext},
    )


from distutils.core import setup
setup( name = 'quicksect'
        , ext_modules=[ Extension("quicksect", sources=["src/quicksect.c"])]
  )
