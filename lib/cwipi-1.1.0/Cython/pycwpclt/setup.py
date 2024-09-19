from distutils.core import setup
from distutils.sysconfig import get_python_lib
import glob
import os
import sys

setup(
    name = "pycwpclt",
    packages     = ['pycwpclt'],
    data_files = [('', ["pycwpclt.so"])],
    author = 'E. Quemerais, B. Andrieu, K. Hoogveld, X. Lamboley',
    description = 'Coupling With Interpolation Parallel Interface - Client',
    license = 'LGPL'
    )

