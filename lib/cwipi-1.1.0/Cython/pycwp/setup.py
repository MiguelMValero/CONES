from distutils.core import setup
from distutils.sysconfig import get_python_lib
import glob
import os
import sys

setup(
    name = "pycwp",
    packages     = ['pycwp'],
    data_files = [('', ["pycwp.so"])],
    author = 'E. Quemerais, B. Andrieu, K. Hoogveld, X. Lamboley',
    description = 'Coupling With Interpolation Parallel Interface',
    license = 'LGPL'
    )

