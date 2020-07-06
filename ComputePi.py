# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 12:29:18 2016

@author: mechd
"""

def ComputePi():
    import math
    import numpy as np
    x = 3
    y = 3
    for i in range(100):
        x = x + math.sin(x)
        pi = float("{0:.12f}".format(x))
        if y != pi:
            y = pi
        else:
            return pi
