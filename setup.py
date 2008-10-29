#from distutils.core import setup
from setuptools import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
#

#print 'remember to cython --convert-range intersection.pyx before  setup.py build_ext --inplace'
import os
os.system('cython --convert-range src/quicksect.pyx')

setup( name = 'quicksect'
        , ext_modules=[ Extension("quicksect", sources=["src/quicksect.pyx"])]
        , cmdclass = {'build_ext': build_ext}
  )

