from distutils.core import setup
from distutils.sysconfig import get_python_lib
import glob
import os
import sys

setup(
    name = "Pypdm",
    packages     = ['Pypdm'],
    data_files = [('', ["Pypdm.so"])],
    author = 'E. Quemerais, B. Maugars, B. Andrieu, K.Hoogveld, J. Coulet',
    author_email='eric.quemerais@onera.fr, bruno.maugars@onera.fr, bastien.andrieu@onera.fr, julien.coulet@onera.fr',
    description = 'Toolkit for parallel distributed computational geometry',
    license = 'LGPL',
    )
