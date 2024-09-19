from distutils.core import setup
from distutils.sysconfig import get_python_lib
import glob
import os
import sys

setup(
    name = "cwipi",
    packages     = ['cwipi'],
    data_files = [('', ["cwipi.so"])],
    author = 'E. Quemerais',
    description = 'Coupling With Interpolation Parallel Interface',
    license = 'LGPL',
    include_dirs=[numpy.get_include()]
    )

