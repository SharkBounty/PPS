from setuptools import setup
from Cython.Build import cythonize
import numpy
import os 

os.chdir('PPS_Python')
# Substitua 'your_module_name' pelo nome que você deu ao seu arquivo .pyx
# O nome do módulo compilado também será 'your_module_name'

setup(
    ext_modules=cythonize("correct_position.pyx"),
    include_dirs=[numpy.get_include()]
)