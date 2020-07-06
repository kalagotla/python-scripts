#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 13:41:36 2018

@author: kalagodk
"""

def diff_flow(gridfile, flowfile1, flowfile2):
    import numpy as np
    import matplotlib.pyplot as plt
    size = 700
    xp, yp = np.zeros(size), np.zeros(size)
    up, vp = np.zeros(size), np.zeros(size)
    g = np.loadtxt(gridfile, skiprows=3)
    f1 = np.loadtxt(flowfile1, skiprows=4)
    f2 = np.loadtxt(flowfile2, skiprows=4)
    f = f1 - f2
    
    kc = 0
    for i in range(size):
        xp[i] = g[kc]
        kc = kc+1
    for i in range(size):
        yp[i] = g[kc]
        kc = kc+1    
    
    kc = 0
    for i in range(size):
        up[i] = f[kc]
        kc = kc+1
    for i in range(size):
        vp[i] = f[kc]
        kc = kc+1
            
    plt.figure(1)
    plt.plot(xp, up, 'g')
    plt.xlabel('X(mm)')
    plt.ylabel('U(m/s)')
    plt.pause(0.05)
    
    plt.figure(2)
    plt.plot(xp, vp, 'g')
    plt.xlabel('X(mm)')
    plt.ylabel('V(m/s)')
    plt.pause(0.05)
    
    return
    