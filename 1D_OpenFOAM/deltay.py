#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 23 15:51:40 2022

@author: miguel
"""

import numpy as np

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

yPlus = 5;
ReTau = 950;
delta = 1;
deltaymin = 9.60e-3;

nu = np.arange(3e-5, 0.0002, 1e-5);

uTau = nu * ReTau / delta;
y = yPlus * nu / uTau;

wished = find_nearest(y, deltaymin);
result = np.where(y == wished);

print('The best value for our initial cell surrounding the Lagrangian markers is : ', float(nu[result]))