# -*- coding: utf-8 -*-
"""
Created on Sat Feb 17 17:15:30 2018

@author: mechd
"""

def wave_equation_g(t, itr):
    import numpy as np
    import matplotlib.pyplot as plt
    import copy
    import scipy.integrate as intg
    from itertools import chain
    
    x1 = np.linspace(-2, -1, itr)
    g1 = np.zeros(itr)
    
    x2 = np.linspace(-1, 0, itr)
    g2 = np.full(itr, -2*t)
    
    x3 = np.linspace(0, 1, itr)
    g3 = np.full(itr, 2*t)
    
    x4 = np.linspace(1, 2, itr)
    g4 = np.zeros(itr)
    
    x_1 = np.append(x1, x2)
    x_2 = np.append(x3, x4)
    x = np.append(x_1, x_2)
    
    g_1 = np.append(g1, g2)
    g_2 = np.append(g3, g4)
    g = np.append(g_1, g_2)
#    
##    print(x)
##    print(g)
    plt.plot(x, g, linewidth=5.0)
    plt.pause(0.05)
    plt.clf()
    plt.ylim(-7, 7)
    
    return