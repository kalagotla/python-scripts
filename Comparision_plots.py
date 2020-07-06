# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 11:57:19 2017

@author: mechd
"""

def postscript_plots(filename, tol, method, psize):
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
    start_time = time.time()
    def simplecount(filename):
        n = 0
        for line in open(filename):
            n += 1
        return n
    n = simplecount(filename)-1
    print('n=', n)
    
    
    def postscript(filename, n, tol):
#    f = open(filename,"r")
#    data = f.readlines()
#    f = np.genfromtxt(filename, delimiter=' ', dtype=None)
        colors = "bgcmyk"
        multi_color = random.choice(colors)
        f = np.loadtxt(filename)
        X1, X2, X3 = np.zeros(n), np.zeros(n), np.zeros(n)
        V1, V2, V3 = np.zeros(n), np.zeros(n), np.zeros(n)
        U1, U2, U3 = np.zeros(n), np.zeros(n), np.zeros(n)
        A1, A2, A3 = np.zeros(n), np.zeros(n), np.zeros(n)
        DT, KC, vj = np.zeros(n), np.zeros(n), np.zeros(n)
        Xf1, Xf2, Xf3 = np.zeros(n), np.zeros(n), np.zeros(n)
        T, axj, ayj = np.zeros(n), np.zeros(n), np.zeros(n)
        X, Y, wj = np.zeros(n), np.zeros(n), np.zeros(n)
        xj, yj, uj = np.zeros(n), np.zeros(n), np.zeros(n)
        """Read data into arrays"""
        for i in range(n):
            X1[i] = f[i, 0]
            X2[i] = f[i, 1]
            X3[i] = f[i, 2]
            V1[i] = f[i, 3]
            V2[i] = f[i, 4]
            V3[i] = f[i, 5]
            U1[i] = f[i, 6]
            U2[i] = f[i, 7]
            U3[i] = f[i, 8]
            A1[i] = f[i, 9]
            A2[i] = f[i, 10]
            A3[i] = f[i, 11]
            DT[i] = f[i, 12]
            KC[i] = f[i, 13]
        """Calculate fluid position"""
        Xf1[0] = X1[0]
        Xf2[0] = X2[0]
        Xf3[0] = X3[0]
        for i in range(n-1):
            Xf1[i+1] = Xf1[i] + DT[i]*U1[i]
            Xf2[i+1] = Xf2[i] + DT[i]*U2[i]
            Xf3[i+1] = Xf3[i] + DT[i]*U3[i]
        """Calculte time for plots"""
        for i in range(n-1):
            T[i+1] = T[i] + DT[i]
        plt.figure(77)
        plt.plot(T, V3)
        return
    return