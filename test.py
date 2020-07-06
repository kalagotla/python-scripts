# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 16:01:01 2016

@author: mechd
"""

def test(v0, theta, g):
    import numpy as np
    #import matplotlib.pyplot as plt
    #N = 100
    #t= np.linspace(0, 1, N)
    #y = np.zeros(N)
    #x = np.zeros(N)
    #for i in range(N-1):
      #  y[i] = v0*np.sin(theta)*t[i] - 0.5*g*t[i]*t[i]
     #   x[i] = v0*np.cos(theta)*t[i]
    #plt.plot(x, y, 'ro')
    
    
    
    for i in range(0, np.pi/2):
        u0 = 10*np.cos(i)
        v0 = 10*np.sin(i)
        i += np.pi/9
    return u0, v0
