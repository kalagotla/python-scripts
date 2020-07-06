# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 14:00:59 2017

@author: mechd
"""
"""
xe, ye, Ue, Ve = experimental_data("M275_theta775_strm.grid.dat", "M275_theta775_strm.req.dat", 35, 86, 1)
"""

def experimental_data(gridfile, flowfile, imax, jmax, kmax):
    import time
    import pdb
    import numpy as np
    import math
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import scipy.interpolate
    from mpl_toolkits.mplot3d import Axes3D
    import random
    from scipy.spatial import Delaunay
    f1 = np.loadtxt(gridfile, skiprows=1)
    f2 = np.loadtxt(flowfile, skiprows=1)
#    x, y = np.zeros((imax, jmax)), np.zeros((imax, jmax))
#    U, V = np.zeros((imax, jmax)), np.zeros((imax, jmax))
    
    x, y = np.zeros(imax*jmax), np.zeros(imax*jmax)
    U, V = np.zeros(imax*jmax), np.zeros(imax*jmax)
    k1, k2 = 0, 0
#    for j in range(jmax):
#        for i in range(imax):
#            if k2>=5:
#                k1 = k1 + 1
#                k2 = 0
#            x[i, j] = f1[k1, k2]
#            U[i, j] = f2[k1, k2]
#            k2 = k2 + 1
#    for j in range(jmax):
#        for i in range(imax):
#            if k2>=5:
#                k1 = k1 + 1
#                k2 = 0
#            y[i, j] = f1[k1, k2]
#            V[i, j] = f2[k1, k2]
#            k2 = k2 + 1
            
    n = 0
    for j in range(jmax):
        for i in range(imax):
            if k2>=5:
                k1 = k1 + 1
                k2 = 0
            x[n] = f1[k1, k2]
            U[n] = f2[k1, k2]
            n = n+1
            k2 = k2 + 1
    n = 0
    for j in range(jmax):
        for i in range(imax):
            if k2>=5:
                k1 = k1 + 1
                k2 = 0
            y[n] = f1[k1, k2]
            V[n] = f2[k1, k2]
            n = n+1
            k2 = k2 + 1            
    
    
    #U = np.fliplr(U)
    #V = np.fliplr(V)
    #print(x, y, U/603.0)
    plt.figure(111)
    size = 700
    
    xi, yi = np.linspace(18.191, 61.851, size), np.linspace(0.25681, 17.721, size)
    xi, yi = np.meshgrid(xi, yi)
        
        # Interpolate; there's also method='cubic' for 2-D data such as here
    zi1 = scipy.interpolate.griddata((x, y), U/603.0, (xi, yi), method='linear')
    
    
    plt.contour(xi, yi, zi1, np.linspace(0, 1.0, 10), linewidths=2.0, linestyles='solid', cmap='jet')
#    plt.imshow(U/603.0, vmin=0, vmax=1.0, origin='lower',
#               extent=[x.min(), x.max(), y.min(), y.max()], cmap='jet')
    #plt.colorbar()
    plt.xlabel('X(mm)')
    plt.ylabel('Y(mm)')
    plt.xlim(18.2, 61.9)
    plt.ylim(0.257, 17.7)
    plt.show()
    plt.figure(112)
    zi2 = scipy.interpolate.griddata((x, y), V/603.0, (xi, yi), method='linear')
    
    plt.contour(xi, yi, zi2, np.linspace(-0.1, 0.1, 10), linewidths=2.0, linestyles='solid', cmap='jet')
#    plt.imshow(V/603.0, vmin=-0.1, vmax=0.1, origin='lower',
#               extent=[x.min(), x.max(), y.min(), y.max()], cmap='jet')
    #plt.colorbar()
    plt.xlabel('X(mm)')
    plt.ylabel('Y(mm)')
    plt.xlim(x.min(), x.max())
    plt.ylim(y.min(), y.max())
    plt.show()
#    return np.fliplr(x), y, U/603.0, V/603.0

    return xi, yi, zi1,zi2