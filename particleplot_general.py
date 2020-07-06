# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 10:32:13 2017

@author: mechd

Can plot normal velocity vs normal distance plot
"""

def particleplotgen(filename, n, tol):
    import pdb
    import numpy as np
    import math
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    import random
#    f = open(filename,"r")
#    data = f.readlines()
#    f = np.genfromtxt(filename, delimiter=' ', dtype=None)
    colors = "bgcmyk"
    multi_color = random.choice(colors)
    f = np.loadtxt(filename, skiprows=1)
    X1, X2, X3 = np.zeros(n), np.zeros(n), np.zeros(n)
    V1, V2, V3 = np.zeros(n), np.zeros(n), np.zeros(n)
    U1, U2, U3 = np.zeros(n), np.zeros(n), np.zeros(n)
    DX1, DX2, DX3 = np.zeros(n), np.zeros(n), np.zeros(n)
    DT, KC, V = np.zeros(n), np.zeros(n), np.zeros(n)
    Xf1, Xf2, Xf3 = np.zeros(n), np.zeros(n), np.zeros(n)
    T, Xt3, Xft3 = np.zeros(n), np.zeros(n), np.zeros(n)
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
    """Creating normalized particle velocity wrt fluid velocity"""
    for i in range(n):
        V[i] = ((round(U3[i], 4) - round(V3[i], 4))
                /(round(U3[i], 4) - round(V3[0], 4)))
        if (np.isnan(V[i]) == True or np.isinf(V[i]) == True
            or V[i] <= 0 or V[i] > 1):
            V[i] = 1.0
    """Moving Shock to zero position"""
    k = 0
    for i in range(n):
        if V[i] < 1.0:
            k = i
            break
#    print(k, X3[k-1])
    if k > 0:
        X3_0 = X3 - X3[k-1]
        """Calculate deflection angle"""
        def_fin = math.atan((X3[n-1]-X3[k-1])/(X1[n-1]-X1[k-1]))
        def_angle = def_fin
        print('theta=', def_angle*180/np.pi)
        deff_fin = math.atan((Xf3[n-1]-Xf3[k-1])/(Xf1[n-1]-Xf1[k-1]))
        deff_angle = deff_fin
        print('thetaf=', deff_angle*180/np.pi)
        def_per = (deff_angle-def_angle)*100/deff_angle
        print('% change in angle = ', def_per)
#    print(len(X1), len(X2), len(X3))
#    Axes3D.plot_wireframe(X1,X2,X3)
#    plt.plot(X1, X2, X3)
#    print(V)
        """Velocity plots"""
    plt.figure(1)
    plt.plot(V1, X2, multi_color, V2, X2, 'r')
    plt.xlabel('U, V(m/s)')
    plt.ylabel('Y(mm)')
    plt.figure(2)
    for i in range(10):
        X2 = X2-i*tol
        Xf2 = Xf2-i*tol
        plt.plot(X1, X2, multi_color, Xf1, Xf2, 'r')
        plt.xlabel('x')
        plt.ylabel('y')
    """Plot xyz in 3-D"""
    fig = plt.figure(3)
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(X1, X2, X3, c='r')
    ax.scatter(Xf1, Xf2, Xf3, c='b')
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    """U, V error"""
    U_error = (sum(V1)-sum(U1))/n
    V_error = (sum(V2)-sum(U2))/n
    print('U_error=', U_error, '\nV_error=', V_error)
    return
    """Plot velocity-contour in X-Y plane"""
    plt.figure()
    X, Y = np.meshgrid(Xf1, Xf2)
    """Can be used to plot avg-vel contour"""
    Vt = np.sqrt(V1**2 + V2**2 + V3**2)
    Vtm, Vx = np.meshgrid(Vt, V1)
#    levels = np.linspace(306,430,100)
#    levels = [306, 308, 310, 312, 320, 350, 370, 400, 405, 415, 420, 425, 430]
#    contour = plt.contour(X, Z, Vz, levels, colors='k') """activate to plot lines with data change"""
#    plt.clabel(contour, colors='k', fmt='%2.1f', fontsize=12)
    contour_filled = plt.pcolor(X, Y, Vtm)
#    plt.contourf(X, Z, Vz, levels) """activate to plot using levels, be careful to chose correct levels"""
    plt.colorbar(contour_filled)
    plt.title('Normal-Velocity Contour in X-Y plane')
    plt.xlabel('x')
    plt.ylabel('Y')
    """Plot fluid-Velocity Contour in X-Y plane"""
    plt.figure()
    Ut = np.sqrt(U1**2 + U2**2 + U3**2)
    Utm, Ux = np.meshgrid(Ut, U1)
#    contour = plt.contour(X, Z, Uz, levels, colors='k', linestyles='dashed')
#    plt.clabel(contour, colors='k', fmt='%2.1f', fontsize=12)
    contour_filled = plt.pcolor(X, Y, Utm)
#    plt.contourf(X, Z, Uz, levels)
    plt.colorbar(contour_filled)
    plt.title('Normal-Fluid-Velocity Contour in X-Y plane')
    plt.xlabel('X')
    plt.ylabel('Y')
    """Plot velocity-contours together to get smeared shock"""
    plt.figure()
#    contourp = plt.contour(X, Z, Vz, levels, colors='k')
#    plt.clabel(contourp, colors='k', fmt='%2.1f', fontsize=12)
    contour_filledp = plt.pcolor(X, Y, Vtm)
#    plt.contourf(X, Z, Vz, levels)
    plt.colorbar(contour_filledp)
    X, Yf = np.meshgrid(Xf1, Xf2-Xf2[n-1]+Xf2[0])
#    contourf = plt.contour(X, Zf, Uz, levels, colors='k', linestyles='dashed')
#    plt.clabel(contourf, colors='k', fmt='%2.1f', fontsize=12)
    contour_filledf = plt.pcolor(X, Yf, Utm)
#    plt.contourf(X, Zf, Uz, levels)
    plt.colorbar(contour_filledf)
    plt.title('Smeared shock Contour in X-Y plane')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.plot(X1, X2, multi_color, Xf1, Xf2-Xf2[n-1]+Xf2[0], 'r')
    plt.show()
    """Plot particle lag vs time"""
    plt.figure()
    for i in range(n-1):
        T[i+1] = T[i] + DT[i]
        Xt3[i+1] = Xt3[i] + V3[i]*DT[i]
        Xft3[i+1] = Xft3[i] + U3[i]*DT[i]
#    print(T[n-1], Xt3-Xft3)
    plt.plot(T, Xt3-Xft3, multi_color)
    plt.title('Particle Lag vs Time')
    plt.xlabel('Time(s)')
    plt.ylabel('Lag behind particle(m)')
    """Plot particle path vs time"""
    plt.figure()
    plt.plot(T, Xt3, multi_color, T, Xft3, multi_color)
    plt.title('Particle path vs Time')
    plt.xlabel('Time(s)')
    plt.ylabel('Distance travelled by particle normal to the shock')
    """Plot slip velocity vs time"""
    plt.figure()
    plt.plot(T, Vt-Ut, multi_color)
    plt.title('Slip Velocity vs Time')
    plt.xlabel('Time(s)')
    plt.ylabel('Slip-Velocity')
    """Plot relative dispersion vs time"""
    plt.figure()
    Xst = np.sqrt((X1-Xf1)**2+(X2-Xf2)**2+(X3-Xf3)**2)
    plt.plot(T, Xst, multi_color)
    plt.title('Relative Dispersion vs Time')
    plt.xlabel('Time(s)')
    plt.ylabel('Relative Dispersion')
    """Streamlines plot"""
    plt.figure()
    Xp, Yp = np.meshgrid(Xf1, Xf2)
    Vx, Vy = np.meshgrid(V1, V2)
    plt.streamplot(Xp, Yp, Vx, Vy, density=0.5, color=Vtm, linewidth=2)
    """Plot Velocity lag vs time"""
    plt.figure()
    Vst = np.sqrt((V1-U1)**2+(V2-U2)**2+(V3-U3)**2)
    Vg, Ug = np.meshgrid(Vst, U3)
    contour_filled = plt.pcolor(X, Y, Vg)
    plt.title('Velocity lag contour')
    plt.xlabel('X')
    plt.ylabel('Y')
#    cp = plt.contour(X, Z, Vz)
#    plt.clabel(cp, inline=True, fontsize=10)
#    var = [float(n) for n in line.split() if n.isdigit() for line in data[1]]
#    print(X1[0], X1[1])
    return