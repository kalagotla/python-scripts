# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 11:05:45 2017

@author: mechd
"""

def contours(filename):
    def simplecount(filename):
        n = 0
        for line in open(filename):
            n += 1
        return n
    n = simplecount(filename)-1
    
    def postscript(filename, n):
        import pdb
        import numpy as np
        import math
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm
        import scipy.interpolate
        from mpl_toolkits.mplot3d import Axes3D
        import random
        colors = "bgcmyk"
        multi_color = random.choice(colors)
        f = np.loadtxt(filename)
        X1, X2, X3 = np.zeros(n), np.zeros(n), np.zeros(n)
        V1, V2, V3 = np.zeros(n), np.zeros(n), np.zeros(n)
        U1, U2, U3 = np.zeros(n), np.zeros(n), np.zeros(n)
        DX1, DX2, DX3 = np.zeros(n), np.zeros(n), np.zeros(n)
        DT, KC, V = np.zeros(n), np.zeros(n), np.zeros(n)
        Xf1, Xf2, Xf3 = np.zeros(n), np.zeros(n), np.zeros(n)
        T, Xt3, Xft3 = np.zeros(n), np.zeros(n), np.zeros(n)
        X, Y = np.zeros(n), np.zeros(n)
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
            DX1[i]= f[i, 9]
            DX2[i]= f[i, 10]
            DX3[i]= f[i, 11]
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
        """Calculate X for plots"""
        X = 25.4*(2.75*X1-38.265)
        Y = 25.4*2.75*X3
        """Use delaunay triangulation for contours"""
        plt.figure(11)
        
        plt.xlim(18, 62)
        plt.ylim(0, 18)
        plt.show()
        return
        
    
    """Call functions for processing"""
    postscript(filename, n, tol)
    return