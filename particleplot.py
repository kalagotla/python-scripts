# -*- coding: utf-8 -*-
"""
Created on Fri May 26 11:02:45 2017

@author: mechd

This script only works for the initial case!!
"""


def particleplot(filename, tol):
#   Relaxation time works correctly for small particles. Use caution for big particles
    import numpy as np
    import math
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    import random
#    f = open(filename,"r")
#    data = f.readlines()
#    f = np.genfromtxt(filename, delimiter=' ', dtype=None)
    n = 0
    for line in open(filename):
        n += 1
    n = n - 1
    print('n=', n)
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
        V[i] = ((round(U3[n-1], 4) - round(V3[i], 4))
                /(round(U3[n-1], 4) - round(V3[0], 4)))
        if (np.isnan(V[i]) == True or np.isinf(V[i]) == True
            or V[i] < 0 or V[i] > 1):
            V[i] = 1.0
        if (V[i] < 1e-5):
            V[i] = -0
    return V, U3, V3
    """Moving Shock to zero position"""
    k = 0
    j = 0
    for i in range(n):
        if V[i] < 1.0:
            k = i
            break
    for i in range(n):
        if V[i] == -0:
            j = i
            break
    print('k=', k)
    print('j=', j)
#    print(k, X3[k-1])
    X3_0 = X3 - X3[k-1]
    """Calculate deflection angle"""
    def_int = math.atan((X3[k-1]-X3[0])/(X1[k-1]-X1[0]))
    def_fin = math.atan((X3[n-1]-X3[k-1])/(X1[n-1]-X1[k-1]))
    def_angle = def_int-def_fin
    print('theta=', def_angle*180/np.pi)
    deff_int = math.atan((Xf3[k-1]-Xf3[0])/(Xf1[k-1]-Xf1[0]))
    deff_fin = math.atan((Xf3[n-1]-Xf3[k-1])/(Xf1[n-1]-Xf1[k-1]))
    deff_angle = deff_int-deff_fin
    print('thetaf=', deff_angle*180/np.pi)
    def_per = (deff_angle-def_angle)*100/deff_angle
    print('% change in angle = ', def_per)
#    print(len(X1), len(X2), len(X3))
#    Axes3D.plot_wireframe(X1,X2,X3)
#    plt.plot(X1, X2, X3)
#    print(V)
    """Velocity plots"""
    plt.figure(1)
    plt.plot(X3_0, V3, multi_color, X3_0, U3, 'r')
    plt.xlabel('Normal distance from shock front(m)')
    plt.ylabel('Normal-Velocity to the shock(m/s)')
    """Calculate time"""
    for i in range(n-1):
        T[i+1] = T[i] + DT[i]
        Xt3[i+1] = Xt3[i] + V3[i]*DT[i]
        Xft3[i+1] = Xft3[i] + U3[i]*DT[i]
    for i in range(10):
        X3 = X3-i*tol
        Xf3 = Xf3-i*tol
        plt.figure(2)
        plt.plot(X1, X3, multi_color, Xf1, Xf3, 'r')
        plt.xlabel('x')
        plt.ylabel('z')
    """Measure particle relaxation time"""
    x = (T[k-1:j-1] - T[k-1])
    y = ((V3[k-1:j-1]-U3[n-1])/(U3[0]-U3[n-1]))
    tau_p = np.polyfit(-x, np.log(y), 1)
    tau_p_linear = sum(DT[k-1:j-1])
    print('tau_p=', 1/tau_p[0])
    print('tau_p_linear=', tau_p_linear)
    plt.figure(3)
    plt.plot(x, y)
    return
    """Plot velocity-contour in X-Y plane"""
    plt.figure(4)
    X, Z = np.meshgrid(Xf1, Xf3)
#    Vt = np.sqrt(V1**2 + V3**2) """Can be used to plot avg-vel contour"""
    Vz, Vx = np.meshgrid(V3, V1)
#    levels = np.linspace(306,430,100)
#    levels = [306, 308, 310, 312, 320, 350, 370, 400, 405, 415, 420, 425, 430]
#    contour = plt.contour(X, Z, Vz, levels, colors='k') """activate to plot lines with data change"""
#    plt.clabel(contour, colors='k', fmt='%2.1f', fontsize=12)
    contour_filled = plt.pcolor(X, Z, Vz)
#    plt.contourf(X, Z, Vz, levels) """activate to plot using levels, be careful to chose correct levels"""
    plt.colorbar(contour_filled)
    plt.title('Normal-Velocity Contour in X-Z plane')
    plt.xlabel('x')
    plt.ylabel('z')
    """Plot fluid-Velocity Contour in X-Y plane"""
    plt.figure(5)
    Uz, Ux = np.meshgrid(U3, U1)
#    contour = plt.contour(X, Z, Uz, levels, colors='k', linestyles='dashed')
#    plt.clabel(contour, colors='k', fmt='%2.1f', fontsize=12)
    contour_filled = plt.pcolor(X, Z, Uz)
#    plt.contourf(X, Z, Uz, levels)
    plt.colorbar(contour_filled)
    plt.title('Normal-Fluid-Velocity Contour in X-Z plane')
    plt.xlabel('x')
    plt.ylabel('z')
    """Plot velocity-contours together to get smeared shock"""
    plt.figure(6)
#    contourp = plt.contour(X, Z, Vz, levels, colors='k')
#    plt.clabel(contourp, colors='k', fmt='%2.1f', fontsize=12)
    contour_filledp = plt.pcolor(X, Z, Vz)
#    plt.contourf(X, Z, Vz, levels)
    plt.colorbar(contour_filledp)
    X, Zf = np.meshgrid(Xf1, Xf3-Xf3[n-1]+Xf3[0])
#    contourf = plt.contour(X, Zf, Uz, levels, colors='k', linestyles='dashed')
#    plt.clabel(contourf, colors='k', fmt='%2.1f', fontsize=12)
    contour_filledf = plt.pcolor(X, Zf, Uz)
#    plt.contourf(X, Zf, Uz, levels)
    plt.colorbar(contour_filledf)
    plt.plot(X1, X3, multi_color, Xf1, Xf3-Xf3[n-1]+Xf3[0], 'r')
    plt.show()
    """Plot particle lag vs time"""
    plt.figure(7)
#    print(T[n-1], Xt3-Xft3)
    plt.plot(T, Xt3-Xft3, multi_color)
    plt.xlabel('Time(s)')
    plt.ylabel('Lag behind particle(m)')
    """Plot particle path vs time"""
    plt.figure(8)
    plt.plot(T, Xt3, multi_color, T, Xft3, multi_color)
    plt.xlabel('Time(s)')
    plt.ylabel('Distance travelled by particle normal to the shock')
    """Plot slip velocity vs time"""
    plt.figure(9)
    plt.plot(T, V3-U3, multi_color)
    plt.xlabel('Time(s)')
    plt.ylabel('Slip-Velocity')
    """Plot relative dispersion vs time"""
    plt.figure(10)
    Xst = np.sqrt((X1-Xf1)**2+(X2-Xf2)**2+(X3-Xf3)**2)
    plt.plot(T, Xst, multi_color)
    plt.xlabel('Time(s)')
    plt.ylabel('Relative Dispersion')
    """Streamlines plot"""
    plt.figure(11)
    plt.streamplot(X, Z, Vx, Vz, color=Uz, linewidth=2)
    """Plot Velocity lag vs time"""
    plt.figure(12)
    Vst = np.sqrt((V1-U1)**2+(V2-U2)**2+(V3-U3)**2)
    Vg, Ug = np.meshgrid(Vst, U3)
    contour_filled = plt.pcolor(X, Z, Vg)
#    cp = plt.contour(X, Z, Vz)
#    plt.clabel(cp, inline=True, fontsize=10)
#    var = [float(n) for n in line.split() if n.isdigit() for line in data[1]]
#    print(X1[0], X1[1])
    return