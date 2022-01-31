from setuptools import setup, find_packages
from Cython.Build import cythonize

setup(name='pybeam',
      version='0.1',
      packages=find_packages(),
      ext_modules = cythonize(["pybeam/default/functions_default.pyx"]),
      zip_safe=False)
