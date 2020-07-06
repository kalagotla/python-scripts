#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 12:20:38 2018

@author: kalagodk
"""

def vel_dia(gridfile, flowfile, i1, j1, dp):
    import numpy as np
    import matplotlib.pyplot as plt
    
    size = 700
    xp, yp = np.zeros((size, size)), np.zeros((size, size))
    up, vp = np.zeros((size, size)), np.zeros((size, size))
    g1 = np.loadtxt(gridfile, skiprows=3)
    f1 = np.loadtxt(flowfile, skiprows=4)
    
    
    kc = 0
    for i in range(size):
        for j in range(size):
            xp[i, j] = g1[kc]
            kc = kc+1
    for i in range(size):
        for j in range(size):
            yp[i, j] = g1[kc]
            kc = kc+1    
    
    kc = 0
    for i in range(size):
        for j in range(size):
            up[i, j] = f1[kc]
            kc = kc+1
    for i in range(size):
        for j in range(size):
            vp[i, j] = f1[kc]
            kc = kc+1
            
    plt.figure(1)
    plt.plot(dp, up[i1, j1], 'o', c='r')
    plt.xlabel('Diameter of Particle(nm)')
    plt.ylabel('Streamwise velocity - Scaled')
    plt.pause(0.05)
    
    plt.figure(2)
    plt.plot(dp, vp[i1, j1], 'o', c='r')
    plt.xlabel('Diameter of Particle(nm)')
    plt.ylabel('Transverse velocity - Scaled')
    plt.pause(0.05)
    
    print('dp=', dp)
    
    return