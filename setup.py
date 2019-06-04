from setuptools.extension import Extension
from setuptools import setup, find_packages

extensions = [Extension("quicksect", ["src/quicksect.pyx"])]

setup(version='0.2.2',
	  name='quicksect',
      description="fast, simple interval intersection",
      long_description=open('README.rst').read(),
      author="Brent Pedersen",
      author_email="bpederse@gmail.com",
      packages=find_packages(),
      setup_requires=['cython'],
      install_requires=['cython'],
      ext_modules=extensions,
      test_suite='nose.collector',
      tests_require='nose',
)
