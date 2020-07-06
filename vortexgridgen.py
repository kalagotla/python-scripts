# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 18:51:39 2018

@author: mechd
"""

def vortexgridgen(n, x, z, r):
    import numpy as np
    import matplotlib.pyplot as plt
#    x = np.linspace(-xmax, xmax, n)
#    y = np.sqrt(r**2 - x**2)
#    y_neg = np.negative(y)
# At zero two points are overlapping
#    y = np.sqrt((r-i)**2 - x**2)
#    y_neg = np.negative(y)
# Plot x, y and y_neg for 2d grid
    
    dx = np.linspace(-x, x, n)
    dy = np.sqrt(r**2 - dx**2)
    
    plt.figure()
    plt.plot(dx, dy, 'ro')
#    bitch please
    
    
    return