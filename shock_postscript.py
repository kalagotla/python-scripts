#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 23:11:17 2019

@author: kalagodk
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 12:12:49 2019

@author: kalagodk
"""

# Code is set for 70% span and 90% pitch
def shock_ps(filename, tol, method, psize, color, logic_f, logic_p):
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
            
        """Calculate zero shock position"""
        for i in range(n):
            if abs(U3[i+1]-U3[i]) >= 0.1:
                shock_loc = Z[i+1]
                break
        """Calculate one particle settling length"""
        for i in range(n):
            if (V3[i]-U3[n]) <= 0.01*U3[n]:
                ps_loc = Z[i]
                break
        print('particle settling length = ', ps_loc)
        
        def plots(color, logic_f, logic_p):
            """Plot streamlines"""
            plt.figure(15)
            if logic_p:
                plt.plot(X, Y, color, label=psize)
                plt.legend()
            if logic_f:
                plt.plot(Xf, Yf, 'b', label='fluid')
                plt.legend()
            #plt.plot(X, V1/603.0, 'g')
            plt.xlabel('X(m)')
            plt.ylabel('Y(m)')
            plt.pause(0.05)
            plt.figure(16)
            if logic_p:
                plt.plot(X, Z, color, label=psize)
                plt.legend()
            if logic_f:
                plt.plot(Xf, Zf, 'b', label='fluid')
                plt.legend()
            #plt.plot(X, V1/603.0, 'g')
            plt.xlabel('X(m)')
            plt.ylabel('Z(m)')
            plt.pause(0.05)
            plt.figure(17)
            if logic_p:
                plt.plot(Y, Z, color, label=psize)
                plt.legend()
            if logic_f:
                plt.plot(Yf, Zf, 'b', label='fluid')
                plt.legend()
            plt.xlabel('Y(m)')
            plt.ylabel('Z(m)')
            plt.pause(0.05)
            #return V2
            
            """Plot xyz in 3-D"""
            fig = plt.figure(18)
#            plt.axes(projection='3d')
            ax = fig.gca(projection='3d')
            if logic_p:
                ax.plot(X, Z, Y, c=color, label=psize)
#                ax.legend()
            if logic_f:
                ax.plot(Xf, Zf, Yf, c='b', label='fluid')
                plt.legend()
            plt.xlabel('X(m)')
            plt.ylabel('Z(m)')
#            plt.zlabel('Y(mm)')
        #plt.zlabel('Z')
            plt.pause(0.05)
        #ax.scatter(Xf1, Xf2, Xf3, c='b')
        
#        return V2
        
            """Velocity plots"""
            plt.figure(19)
            if logic_p:
                plt.plot(Z-shock_loc, V1, color, label=psize)
                plt.legend()
            if logic_f:
                plt.plot(Z-shock_loc, U1, 'b', label='fluid')
                plt.legend()
            plt.xlabel('X(m)')
            plt.ylabel('U(m/s)')
            plt.pause(0.05)
            plt.figure(20)
            if logic_p:
                plt.plot(Z-shock_loc, V3, color, label=psize)
                plt.legend()
            if logic_f:
                plt.plot(Z-shock_loc, U3, 'b', label='fluid')
                plt.legend()
            plt.xlabel('Shock Normal Distance(m)')
            plt.ylabel('Shock Normal Velocity(m/s)')
            plt.pause(0.05)
            plt.figure(21)
            if logic_p:
                plt.plot(Z-shock_loc, V2, color, label=psize)
                plt.legend()
            if logic_f:
                plt.plot(Z-shock_loc, U2, 'b', label='fluid')
                plt.legend()
            plt.xlabel('X(m)')
            plt.ylabel('W(m/s)')
            plt.pause(0.05)
            """In terms of R37 chord length"""
            plt.figure(22)
            if logic_p:
                plt.plot((Z-shock_loc)*100/0.055905, V3, color, label=psize)
                plt.legend()
            if logic_f:
                plt.plot((Z-shock_loc)*100/0.055905, U3, 'b', label='fluid')
                plt.legend()
            plt.xlabel('% Chord Length')
            plt.ylabel('Shock Normal Velocity(m/s)')
            plt.pause(0.05)
            print('PSL in terms of CL = ', (ps_loc-shock_loc)*100/0.055905)
            
            return
        
        plots(color, logic_f, logic_p)
        
        
        return
    
    
    single_line_plots(filename, n, tol, color, logic_f, logic_p)
    end_time = time.time()
    print('Elapsed_time = ', end_time-start_time)
    return