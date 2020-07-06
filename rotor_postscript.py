#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 12:12:49 2019

@author: kalagodk
"""

# Code is set for 70% span and 90% pitch
def rotor_ps(filename, tol, method, psize, color, logic_f, logic_p):
    import time
    import pdb
    import numpy as np
    import math
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import scipy.interpolate
    from mpl_toolkits.mplot3d import Axes3D
    import random
    import Experimental_data as ed
    import CFD as cfd
    from scipy.spatial import Delaunay
    start_time = time.time()
    def simplecount(filename):
        n = 0
        for line in open(filename):
            n += 1
        return n
    n = simplecount(filename)-1
    print('n=', n)
    
    def single_line_plots(filename, n, tol, color, logic_f, logic_p):
        f = np.loadtxt(filename)
        X, Y, Z       = np.zeros(n), np.zeros(n), np.zeros(n)
        V1, V2, V3    = np.zeros(n), np.zeros(n), np.zeros(n)
        U1, U2, U3    = np.zeros(n), np.zeros(n), np.zeros(n)
        A1, A2, A3    = np.zeros(n), np.zeros(n), np.zeros(n)
        Av1, Av2, Av3 = np.zeros(n), np.zeros(n), np.zeros(n)
        Au1, Au2, Au3 = np.zeros(n), np.zeros(n), np.zeros(n)
        DT, KC        = np.zeros(n), np.zeros(n)
        Xf, Yf, Zf    = np.zeros(n), np.zeros(n), np.zeros(n)
        """Load respective variables from file"""
        X  = f[:, 0]
        Y  = f[:, 1]
        Z  = f[:, 2]
        V1 = f[:, 3]
        V2 = f[:, 4]
        V3 = f[:, 5]
        U1 = f[:, 6]
        U2 = f[:, 7]
        U3 = f[:, 8]
#        A1 = f[:, 9]
#        A2 = f[:,10]
#        A3 = f[:,11]
#        Av1= f[:,12]
#        Av2= f[:,13]
#        Av3= f[:,14]
#        Au1= f[:,15]
#        Au2= f[:,16]
#        Au3= f[:,17]
#        DT = f[:,18]
#        KC = f[:,19]
        """Calculate fluid position"""
        Xf[0] = X[0]
        Yf[0] = Y[0]
        Zf[0] = Z[0]
        for i in range(n-1):
            Xf[i+1] = Xf[i] + DT[i]*U1[i]
            Yf[i+1] = Yf[i] + DT[i]*U2[i]
            Zf[i+1] = Zf[i] + DT[i]*U3[i]
        """Calculate % Chord"""
        cp, temp = np.zeros(n), np.zeros(n)
        cl = 0.055905 # Chord length
        ymin, ymax = -0.02385545, 0.0220293
        zmin, zmax =  0.00427579, 0.0362121
        theta = np.arcsin((ymax - ymin)/cl)
        l = 0.02454 # Rotor 37 pitch
        # End coordinates of chord at 90% pitch and 70% span
        y10max = ymax - 0.1*l*np.sin(theta)
        y10min = ymin - 0.1*l*np.sin(theta)
        z10max = np.sqrt((0.1*l)**2 - (y10max-ymax)**2)
        z10min = np.sqrt((0.1*l)**2 - (y10min-ymin)**2)
        temp = np.sign(Y-y10min)
        cp = temp*np.sqrt((Y-y10min)**2+(Z-z10min)**2)*100/cl
        
        """Velocity Magnitudes"""
        Vm, Um = np.zeros(n), np.zeros(n)
        Vm = np.sqrt(V1**2+V2**2+V3**2)
        Um = np.sqrt(U1**2+U2**2+U3**2)
        
        def plots(color, logic_f, logic_p):
            """Plot streamlines"""
            plt.figure(15)
            if logic_p:
                plt.plot(X, Y, color, label=psize)
                plt.legend()
            if logic_f:
                plt.plot(Xf, Yf, 'b', label='liquid')
                plt.legend()
            plt.xlabel('X(m)')
            plt.ylabel('Y(m)')
            plt.pause(0.05)
            plt.figure(16)
            if logic_p:
                plt.plot(X, Z, color, label=psize)
                plt.legend()
            if logic_f:
                plt.plot(Xf, Zf, 'b', label='liquid')
                plt.legend()
            plt.xlabel('X(m)')
            plt.ylabel('Z(m)')
            plt.pause(0.05)
            plt.figure(17)
            if logic_p:
                plt.plot(Y, Z, color, label=psize)
                plt.legend()
            if logic_f:
                plt.plot(Yf, Zf, 'b', label='liquid')
                plt.legend()
            plt.xlabel('Y(m)')
            plt.ylabel('Z(m)')
            plt.pause(0.05)
            
            """Plot xyz in 3-D"""
            fig = plt.figure(18)
            ax = fig.gca(projection='3d')
            if logic_p:
                ax.plot(X, Z, Y, c=color, label=psize)
            if logic_f:
                ax.plot(Xf, Zf, Yf, c='b', label='liquid')
                plt.legend()
            plt.xlabel('X(m)')
            plt.ylabel('Z(m)')
            plt.pause(0.05)
        
            """Velocity plots"""
            plt.figure(19)
            if logic_p:
                plt.plot(cp, V1, color, label=psize)
                plt.legend()
            if logic_f:
                plt.plot(cp, U1, 'b', label='liquid')
                plt.legend()
            plt.xlabel('% Chord')
            plt.ylabel('U(m/s)')
            plt.pause(0.05)
            plt.figure(20)
            if logic_p:
                plt.plot(cp, V3, color, label=psize)
                plt.legend()
            if logic_f:
                plt.plot(cp, U3, 'b', label='liquid')
                plt.legend()
            plt.xlabel('% chord')
            plt.ylabel('W(m/s)')
            plt.pause(0.05)
            plt.figure(21)
            if logic_p:
                plt.plot(cp, V2, color, label=psize)
                plt.legend()
            if logic_f:
                plt.plot(cp, U2, 'b', label='liquid')
                plt.legend()
            plt.xlabel('% chord')
            plt.ylabel('V(m/s)')
            plt.pause(0.05)
            plt.figure(22)
            if logic_p:
                plt.plot(cp, Vm, color, label=psize)
                plt.legend()
            if logic_f:
                plt.plot(cp, Um, 'b', label='liquid')
                plt.legend()
            plt.xlabel('% chord')
            plt.ylabel('Velocity magnitude(m/s)')
            plt.pause(0.05)
            
            return
        
        plots(color, logic_f, logic_p)
        
        
        return
    
    
    single_line_plots(filename, n, tol, color, logic_f, logic_p)
    end_time = time.time()
    print('Elapsed_time = ', end_time-start_time)
    return