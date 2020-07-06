# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 13:09:10 2017

@author: mechd
"""
# xc, yc, Uc, Vc = CFD("CFD_Grid_z0.txt", "CFD_Flow_z0.txt", 1e-4)
def CFD(gridfile, flowfile, tol):
    import time
    import pdb
    import numpy as np
    import math
    import matplotlib.pyplot as plt
    import scipy.interpolate
    from mpl_toolkits.mplot3d import Axes3D
    import random
    import matplotlib.tri as tri
    import pylab
    start_time = time.time()
    def simplecount(gridfile):
        n = 0
        for line in open(gridfile):
            n += 1
        return n
    n = simplecount(gridfile)-1
    print('n=', n)
    
    def post_script(gridfile, flowfile):
        f1 = np.loadtxt(gridfile)
        f2 = np.loadtxt(flowfile)
        X1, X2, X3 = np.zeros(n), np.zeros(n), np.zeros(n)
        U1, U2, U3 = np.zeros(n), np.zeros(n), np.zeros(n)
        X, Y = np.zeros(n), np.zeros(n)
        xj, yj, uj = np.zeros(n), np.zeros(n), np.zeros(n)
        vj, wj = np.zeros(n), np.zeros(n)
        for i in range(n):
            X1[i] = f1[i, 0]
            X2[i] = f1[i, 1]
            X3[i] = f1[i, 2]
            U1[i] = f2[i, 0]
            U2[i] = f2[i, 1]
            U3[i] = f2[i, 2]
        """Calculate X for plots"""
        X = 25.4*(2.75*X1-38.265)
        Y = 25.4*2.75*X3
        Z = 25.4*(2.75*X2-1.125)
        """Isolate experimental plane data"""
        j = 0
        for i in range(n):
            if 17<=X[i]<=65:
                if 0.01<=Y[i]<=20:
                    xj[j] = X[i]
                    yj[j] = Y[i]
                    uj[j] = U1[i]
                    vj[j] = U3[i]
                    wj[j] = U2[i]
                    j=j+1
        print('j=', j)
        xj = np.trim_zeros(xj,'b')
        yj = np.trim_zeros(yj,'b')
        vj = np.trim_zeros(vj,'b')
        uj = np.trim_zeros(uj,'b')
        wj = np.trim_zeros(wj,'b')
        size = 700
        #print('sizes=', np.size(xj), np.size(yj), np.size(uj), np.size(vj))
        """Contour for U-velocity"""
        plt.figure(11)
        #xi, yi = np.linspace(xj.min(), xj.max(), size), np.linspace(yj.min(), yj.max(), size)
        xi, yi = np.linspace(18.191, 61.851, size), np.linspace(0.25681, 17.721, size)
        xi, yi = np.meshgrid(xi, yi)
        
        # Interpolate; there's also method='cubic' for 2-D data such as here
        zi1 = scipy.interpolate.griddata((xj, yj), uj/603.0, (xi, yi), method='linear')
        #print(xi, yi, zi)        
        plt.contour(xi ,yi, zi1, np.linspace(0, 1.0, 10), linewidths=2.0, linestyles='solid', cmap='jet')        
        #plt.imshow(zi1, vmin=0, vmax=1.0, origin='lower',
         #          extent=[xj.min(), xj.max(), yj.min(), yj.max()])
        #triang = tri.Triangulation(xj, yj)
        #plt.tricontourf(triang, uj/603.0)
        #plt.colorbar()
        #plt.xlabel('X(mm)')
        #plt.ylabel('Y(mm)')
        #plt.xlim(18.2, 61.9)
        #plt.ylim(0.257, 17.7)
        plt.title('U-velocity contour')
        #plt.legend()
        """Contour for V-velocity"""
        plt.figure(12)

        # Interpolate; there's also method='cubic' for 2-D data such as here
        zi2 = scipy.interpolate.griddata((xj, yj), vj/603.0, (xi, yi), method='linear')
        
        plt.contour(xi ,yi, zi2, np.linspace(-0.1, 0.1, 10), linewidths=2.0, linestyles='solid', cmap='jet')
        #plt.imshow(zi2, vmin=-0.1, vmax=0.1, origin='lower',
         #      extent=[xj.min(), xj.max(), yj.min(), yj.max()])
        #plt.colorbar()
        #plt.xlabel('X(mm)')
        #plt.ylabel('Y(mm)')
        #plt.xlim(18.2, 61.9)
        #plt.ylim(0.257, 17.7)
        plt.title('V-velocity contour')
        #plt.legend()
        plt.show()
        
        """W-velocity"""
        zi3 = scipy.interpolate.griddata((xj, yj), wj/603.0, (xi, yi), method='linear')
        
        """Write files for tecplot (Plot3D format)"""
        """kc = 0
        #first number: dimensions, 2nd represents imax, jmax, kmax
        #first number: No. of variables, 2nd: imax, jmax, kmax, Max_Variable_no.
        f = np.zeros(size*size*3+3)
        g = np.zeros(size*size*3+4)
        f[kc] = int(size)
        g[kc] = int(size)
        f[kc+1] = int(size)
        g[kc+1]  = int(size)
        f[kc+2] = 1
        g[kc+2] = 1
        g[kc+3] = 3
        kc = 3
        for i in range(size):
            for j in range(size):
                f[kc] = xi[i, j]
                kc = kc+1
        for i in range(size):
            for j in range(size):
                f[kc] = yi[i, j]
                kc = kc+1
        for i in range(size):
            for j in range(size):
                f[kc] = 0
                kc = kc+1
        kc = 4
        for m in range(3):
            for i in range(size):
                for j in range(size):
                    if m==0:
                        g[kc] = zi1[i, j]
                    if m==1:
                        g[kc] = zi2[i, j]
                    if m==2:
                        g[kc] = zi3[i, j]
                    kc = kc+1
            #print(f, g)
        np.savetxt('grid_data_CFD.txt', f)
        np.savetxt('flow_data_CFD.txt', g)"""
        
        
        return xi, yi, zi1, zi2
    
    xc, yc, Uc, Vc = post_script(gridfile, flowfile)
    end_time = time.time()
    print('Elapsed_time=', end_time-start_time)
    return xc, yc, Uc, Vc