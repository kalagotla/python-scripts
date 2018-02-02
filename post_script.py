# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 22:36:55 2017

@author: mechd
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 10:32:13 2017

@author: mechd

# postscript_main("merged_data.txt", 1e-4)
#import matplotlib.pyplot as plt
#plt.figure(16)
#plt.axes(projection='3d')
#for i in range(98):
#    postscript_main("postdata\Z" + str(i+1) + ".rtf", 1e-4)

#xi, yi, ui, vi, wi = postscript_main('../Statistical_distribution/254/merged_data_254.txt', 1e-4, 'linear', 254)
"""

def postscript_main(filename, tol, method, psize):
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
        Av1, Av2, Av3 = np.zeros(n), np.zeros(n), np.zeros(n)
        Au1, Au2, Au3 = np.zeros(n), np.zeros(n), np.zeros(n)
        axvj, ayvj, azvj = np.zeros(n), np.zeros(n), np.zeros(n)
        axuj, ayuj, azuj = np.zeros(n), np.zeros(n), np.zeros(n)
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
            Av1[i] = f[i, 12]
            Av2[i] = f[i, 13]
            Av3[i] = f[i, 14]
            Au1[i] = f[i, 15]
            Au2[i] = f[i, 16]
            Au3[i] = f[i, 17]
            DT[i] = f[i, 18]
            KC[i] = f[i, 19]
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

        X = 25.4*(2.75*X1*1e3-38.265)
        Y = 25.4*2.75*X3*1e3
        Z = 25.4*(2.75*X2*1e3-1.125)
        Xf = 25.4*(2.75*Xf1*1e3-38.265)
        Yf = 25.4*2.75*Xf3*1e3
        Zf = 25.4*(2.75*Xf2*1e3-1.125)
        
        """Acceleration from velocity derivatives"""
        #Cannot do this because you need to do this for each particle separately
        #It's better to do this in Visual3 Streamlines.f scripts
#        ax1 = [(t - s)/1e-9 for s, t in zip(V1, V1[1:])]
#        ay1 = [(t - s)/1e-9 for s, t in zip(V3, V3[1:])]
#        ax1 = np.append(ax1, A1[n-1])
#        ay1 = np.append(ay1, A3[n-1]) 
        #X, Y, Z = X1*1e3, X3*1e3, X2*1e3
        #Xf, Yf, Zf = Xf1*1e3, Xf3*1e3, Xf2*1e3
            
            
        """Use this with individual particles"""
        def plots():
            """Plot streamlines"""
            plt.figure(15)
            plt.plot(X, Y, 'r')
            plt.plot(Xf, Yf, 'b')
            #plt.plot(X, V1/603.0, 'g')
            plt.xlabel('X')
            plt.ylabel('Y')
            plt.pause(0.05)
            plt.figure(16)
            plt.plot(X, Z, 'r')
            plt.plot(Xf, Zf, 'b')
            #plt.plot(X, V1/603.0, 'g')
            plt.xlabel('X')
            plt.ylabel('Z')
            plt.pause(0.05)
            plt.figure(17)
            plt.plot(Y, Z, 'r')
            plt.plot(Yf, Zf, 'b')
            plt.xlabel('Y')
            plt.ylabel('Z')
            plt.pause(0.05)
            #return V2
            
            """Plot xyz in 3-D"""
            plt.figure(18)
            plt.axes(projection='3d')
            plt.plot(X, Z, Y, c='r')
            plt.plot(Xf, Zf, Yf, c='b')
            plt.xlabel('X')
            plt.ylabel('Z')
        #plt.zlabel('Z')
            plt.pause(0.05)
        #ax.scatter(Xf1, Xf2, Xf3, c='b')
        
#        return V2
        
            """Velocity plots"""
            plt.figure(19)
            plt.plot(X, V1, 'r', X, U1, 'b')
            plt.xlabel('X(mm)')
            plt.ylabel('U(m/s)')
            plt.pause(0.05)
            plt.figure(20)
            plt.plot(X, V3, 'r', X, U3, 'b')
            plt.xlabel('X(mm)')
            plt.ylabel('V(m/s)')
            plt.pause(0.05)
            plt.figure(21)
            plt.plot(X, V2, 'r', X, U2, 'b')
            plt.xlabel('X(mm)')
            plt.ylabel('W(m/s)')
            plt.pause(0.05)
            
        
            return
        
        
        
        """Isolate experimental data plane"""
        j = 0
        for i in range(n):
             if 17.7<=X[i]<=65:
                 if 0.01<=Y[i]<=20:
                     xj[j] = X[i]
                     yj[j] = Y[i]
                     uj[j] = V1[i]
                     vj[j] = V3[i]
                     wj[j] = V2[i]
                     axj[j] = A1[i]
                     ayj[j] = A3[i]
                     axvj[j] = Av1[i]
                     ayvj[j] = Av3[i]
                     azvj[j] = Av2[i]
                     axuj[j] = Au1[i]
                     ayuj[j] = Au3[i]
                     azuj[j] = Au2[i]
                     j=j+1
        print('j=', j)
        xj = np.trim_zeros(xj)
        yj = np.trim_zeros(yj)
        vj = np.trim_zeros(vj)
        uj = np.trim_zeros(uj)
        wj = np.trim_zeros(wj)
        axj = np.trim_zeros(axj)
        ayj = np.trim_zeros(ayj)
        axvj = np.trim_zeros(axvj)
        ayvj = np.trim_zeros(ayvj)
        azvj = np.trim_zeros(azvj)
        axuj = np.trim_zeros(axuj)
        ayuj = np.trim_zeros(ayuj)
        azuj = np.trim_zeros(azuj)
        
        
        #plt.figure(15)
        #plt.plot(axj, ax1)
        #plt.xlabel('axj')
        #plt.ylabel('ax1')
        #plt.pause(0.05)
        
        #return
        
        #return xj, yj, uj, vj, wj
        
        #plt.tripcolor(xj, yj, uj)
        #plt.tricontourf(xj, yj, uj, 100)
        #return
        
        def contours():
            size = 700
            """Plot all contours on same plot, U-velocity"""
            plt.figure(11)
            xi, yi = np.linspace(18.2, 61.9, size), np.linspace(0.257, 17.7, size)
            xi, yi = np.meshgrid(xi, yi)
        
        # Interpolate; there's also method='cubic' for 2-D data such as here
            zi1 = scipy.interpolate.griddata((xj, yj), uj/603.0, (xi, yi), method=method, fill_value=1.1)
        #print(xi, yi, zi)        
            plt.contour(xi ,yi, zi1, np.linspace(0, 1.0, 10), linewidths=1.0, linestyles='dashed', cmap='jet')        
            #plt.imshow(zi1, vmin=0, vmax=1, origin='lower',
             #           extent=[xj.min(), xj.max(), yj.min(), yj.max()], cmap='jet')
            plt.colorbar()
            plt.xlabel('X(mm)')
            plt.ylabel('Y(mm)')
            plt.xlim(18.2, 61.9)
            plt.ylim(0.257, 17.7)
            plt.pause(0.05)
        #plt.xlim(xj.min(), xj.max())
        #plt.ylim(yj.min(), yj.max())
            plt.title('U-velocity contour')
            """Contour for V-velocity"""
            plt.figure(12)

        # Interpolate; there's also method='cubic' for 2-D data such as here
            zi2 = scipy.interpolate.griddata((xj, yj), vj/603.0, (xi, yi), method=method, fill_value=1.1)
        
            plt.contour(xi ,yi, zi2, np.linspace(-0.1, 0.1, 10), linewidths=1.0, linestyles='dashed', cmap='jet')
            #plt.imshow(zi2, vmin=-0.1, vmax=0.1, origin='lower',
             #          extent=[xj.min(), xj.max(), yj.min(), yj.max()], cmap='jet')
            plt.colorbar()
            plt.xlabel('X(mm)')
            plt.ylabel('Y(mm)')
            plt.xlim(18.2, 61.9)
            plt.ylim(0.257, 17.7)
            plt.pause(0.05)
        #plt.xlim(xj.min(), xj.max())
        #plt.ylim(yj.min(), yj.max())
            plt.title('V-velocity contour')
            plt.show()
            
            # compare with experimental data
            plt.figure(111)
            plt.contour(xi ,yi, zi1, np.linspace(0, 1.0, 10), linewidths=1.0, linestyles='dashed', cmap='jet')
            plt.colorbar()
            plt.xlabel('X(mm)')
            plt.ylabel('Y(mm)')
            plt.xlim(18.2, 61.9)
            plt.ylim(0.257, 17.7)
            plt.pause(0.05)
            plt.title('U-velocity contour')
            plt.show()
            
            plt.figure(112)
            plt.contour(xi ,yi, zi2, np.linspace(-0.1, 0.1, 10), linewidths=1.0, linestyles='dashed', cmap='jet')
            plt.colorbar()
            plt.xlabel('X(mm)')
            plt.ylabel('Y(mm)')
            plt.xlim(18.2, 61.9)
            plt.ylim(0.257, 17.7)
            plt.pause(0.05)
            plt.title('V-velocity contour')
            plt.show()
        # Z-direction velocity
            zi3 = scipy.interpolate.griddata((xj, yj), wj/603.0, (xi, yi), method=method, fill_value=1.1)
            
            """Write files for tecplot (Plot3D format)"""
            kc = 0
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
            np.savetxt('grid_data_'+ str(psize) +'.txt', f)
            np.savetxt('flow_data_'+ str(psize) +'.txt', g)
                    
        
            #return xi, yi, zi1, zi2, zi3
            """X-acceleration contour"""
            plt.figure(13)

        # Interpolate; there's also method='cubic' for 2-D data such as here
            zi = scipy.interpolate.griddata((xj, yj), axj, (xi, yi), method=method, fill_value=1.1)

            plt.imshow(zi, vmin=axuj.min(), vmax=axuj.max(), origin='lower',
                       extent=[xj.min(), xj.max(), yj.min(), yj.max()], cmap='jet')
            plt.colorbar()
            plt.xlabel('X(mm)')
            plt.ylabel('Y(mm)')
            plt.pause(0.05)
        #plt.xlim(18.2, 61.9)
        #plt.ylim(0.257, 17.7)
            plt.show()
            """Y-acceleration contour"""
            plt.figure(14)

        # Interpolate; there's also method='cubic' for 2-D data such as here
            zi = scipy.interpolate.griddata((xj, yj), ayj, (xi, yi), method=method, fill_value=1.1)
            
            plt.imshow(zi, vmin=ayj.min(), vmax=ayj.max(), origin='lower',
                       extent=[xj.min(), xj.max(), yj.min(), yj.max()], cmap='jet')
            plt.colorbar()
            plt.xlabel('X(mm)')
            plt.ylabel('Y(mm)')
            plt.pause(0.05)
        #plt.xlim(18.2, 61.9)
        #plt.ylim(0.257, 17.7)
            plt.show()
            
            """X-acceleration contour"""
            plt.figure(131)

        # Interpolate; there's also method='cubic' for 2-D data such as here
            zi = scipy.interpolate.griddata((xj, yj), axvj, (xi, yi), method=method, fill_value=1.1)

            plt.imshow(zi, vmin=axuj.min(), vmax=axuj.max(), origin='lower',
                       extent=[xj.min(), xj.max(), yj.min(), yj.max()], cmap='jet')
            plt.colorbar()
            plt.xlabel('X(mm)')
            plt.ylabel('Y(mm)')
            plt.pause(0.05)
        #plt.xlim(18.2, 61.9)
        #plt.ylim(0.257, 17.7)
            plt.show()
            """Y-acceleration contour"""
            plt.figure(141)

        # Interpolate; there's also method='cubic' for 2-D data such as here
            zi = scipy.interpolate.griddata((xj, yj), ayvj, (xi, yi), method=method, fill_value=1.1)
            
            plt.imshow(zi, vmin=ayj.min(), vmax=ayj.max(), origin='lower',
                       extent=[xj.min(), xj.max(), yj.min(), yj.max()], cmap='jet')
            plt.colorbar()
            plt.xlabel('X(mm)')
            plt.ylabel('Y(mm)')
            plt.pause(0.05)
        #plt.xlim(18.2, 61.9)
        #plt.ylim(0.257, 17.7)
            plt.show()
            
            """X-acceleration contour"""
            plt.figure(132)

        # Interpolate; there's also method='cubic' for 2-D data such as here
            zi = scipy.interpolate.griddata((xj, yj), axuj, (xi, yi), method=method, fill_value=1.1)

            plt.imshow(zi, vmin=axuj.min(), vmax=axuj.max(), origin='lower',
                       extent=[xj.min(), xj.max(), yj.min(), yj.max()], cmap='jet')
            plt.colorbar()
            plt.xlabel('X(mm)')
            plt.ylabel('Y(mm)')
            plt.pause(0.05)
        #plt.xlim(18.2, 61.9)
        #plt.ylim(0.257, 17.7)
            plt.show()
            """Y-acceleration contour"""
            plt.figure(142)

        # Interpolate; there's also method='cubic' for 2-D data such as here
            zi = scipy.interpolate.griddata((xj, yj), ayuj, (xi, yi), method=method, fill_value=1.1)
            
            plt.imshow(zi, vmin=ayj.min(), vmax=ayj.max(), origin='lower',
                       extent=[xj.min(), xj.max(), yj.min(), yj.max()], cmap='jet')
            plt.colorbar()
            plt.xlabel('X(mm)')
            plt.ylabel('Y(mm)')
            plt.pause(0.05)
        #plt.xlim(18.2, 61.9)
        #plt.ylim(0.257, 17.7)
            plt.show()
            
            return xi, yi, zi1, zi2, zi3
        
        """Call plots, contours"""
#        plots()
        xi, yi, zi1, zi2, zi3 = contours()
        return xi, yi, zi1, zi2, zi3
        #return
        
        
        
        """-----------Contours above, plots below--------------"""
        
        
        """Plot X-Y for analysis"""
        plt.figure(1)
        plt.plot(X, Y)
        #plt.scatter(X, Y, c=V1/603.0)
        plt.xlim(18.2, 61.9)
        plt.ylim(0.257, 17.7)
        """Velocity plots"""
        plt.figure(2)
        plt.plot(X, V1, multi_color, X, U1, 'r')
        plt.xlabel('X(mm)')
        plt.ylabel('U(m/s)')
        plt.figure(3)
        plt.scatter(X, V2, c=V2, cmap=cm.jet)
        plt.xlabel('X(mm)')
        plt.ylabel('V(m/s)')
        plt.figure(4)
        plt.plot(X, V3, multi_color, X, U3, 'r')
        plt.xlabel('X(mm)')
        plt.ylabel('W(m/s)')
        """Plot particle-fluid together"""
        plt.figure(5)
        for i in range(10):
            X2 = X2-i*tol
            Xf2 = Xf2-i*tol
            plt.plot(X1, X2, multi_color, Xf1, Xf2, 'r')
            plt.xlabel('x')
            plt.ylabel('y')
        """Plot relative dispersion vs X-position"""
        plt.figure(6)
        Xst = np.sqrt((X1-Xf1)**2+(X2-Xf2)**2+(X3-Xf3)**2)
        plt.plot(X, Xst, multi_color)
        plt.title('Relative Dispersion vs Time')
        plt.xlabel('X(mm)')
        plt.ylabel('Relative Dispersion')
        """Plot slip velocity vs X-position"""
        plt.figure(7)
        Vt = np.sqrt(V1**2 + V2**2 + V3**2)
        Ut = np.sqrt(U1**2 + U2**2 + U3**2)
        plt.plot(X, Vt-Ut, multi_color)
        plt.title('Slip Velocity vs Time')
        plt.xlabel('X(mm)')
        plt.ylabel('Slip-Velocity(m/s)')
        """Particle slip in each direction"""
        plt.figure(8)
        plt.plot(T, V1-U1, multi_color)
        plt.title('Slip Velocity X vs Time')
        plt.xlabel('Time(s)')
        plt.ylabel('Slip-Velocity(m/s)')
        plt.figure(9)
        plt.plot(T, V2-U2, multi_color)
        plt.title('Slip Velocity Y vs Time')
        plt.xlabel('Time(s)')
        plt.ylabel('Slip-Velocity(m/s)')
        plt.figure(10)
        plt.plot(T, V3-U3, multi_color)
        plt.title('Slip Velocity Z vs Time')
        plt.xlabel('Time(s)')
        plt.ylabel('Slip-Velocity(m/s)')
        return
    
        """U, V error"""
        U_error = (sum(V1)-sum(U1))/n
        V_error = (sum(V2)-sum(U2))/n
        print('U_error=', U_error, '\nV_error=', V_error)
        
    """Call functions for processing"""
    xi, yi, zi1, zi2, zi3 = postscript(filename, n, tol)
    #x, y, u, v, w = postscript(filename, n, tol)
    #V1 = postscript(filename, n, tol)
    #postscript(filename, n, tol)
    end_time = time.time()
    print('Elapsed_time=', end_time-start_time)
    #return
    return xi, yi, zi1, zi2, zi3
    #return x, y, u, v, w
#    cp = plt.contour(X, Z, Vz)
#    plt.clabel(cp, inline=True, fontsize=10)
#    var = [float(n) for n in line.split() if n.isdigit() for line in data[1]]
#    print(X1[0], X1[1])