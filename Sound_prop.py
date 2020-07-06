# -*- coding: utf-8 -*-
"""
Created on Mon Apr  2 08:14:54 2018

@author: mechd
"""

def sound_prop(f, amp, t_tol, itr, n, m):
    import numpy as np
    import matplotlib.pyplot as plt
    import scipy.interpolate
    import random
    t = np.linspace(0, t_tol, itr)
    x = np.linspace(-m, m, itr)
    y = np.linspace(-m, m, itr)
    c = 340 #Speed of Sound
    k = 2*np.pi*f/c #Wave Number
    omega = 2*np.pi*f #Angular frequency
    rho = 1.225 #Air density
    xs, ys = np.zeros(n), np.zeros(n)
    xs[0] = (x.min()+x.max())/2
    ys[0] = (y.min()+y.max())/2
    i = 1
    for i in range(n):
        xs[i] = random.choice(x)
        ys[i] = random.choice(y)
#        xs[i] = xs[i-1]+(x.max()-x.min())/itr
#        ys[i] = ys[i-1]+(y.max()-y.min())/itr
#    print(xs)
    r = np.zeros(n)
    xi, yi = np.meshgrid(x, y)
    P1 = np.zeros((n, itr, itr))
    P2 = np.zeros((itr, itr))
    plt.figure(1)
    for i1 in t:
        i = 0
        for i2 in x:
            j = 0
            for i3 in y:
                for i4 in range(n):
                    r[i4] = ((i2-xs[i4])**2+(i3-ys[i4])**2)**0.5
                    if r[i4]==0:
                        P1[i4,i,j] = 0.01*amp*np.sin(omega*i1)
#                        print(i4)
                    elif i4<0.5*n:
                        P1[i4,i,j] = (rho/(4*np.pi*r[i4]))*amp*(np.sin(omega*(i1-r[i4]/c)-k*r[i4]))
#                        print('yes1', i4, 0.5*n)
                    else:
                        P1[i4,i,j] = (rho/(4*np.pi*r[i4]))*amp*(np.sin(omega*(i1-r[i4]/c)-k*r[i4]+np.pi))
#                        print('yes2', i4, 0.5*n)
                j += 1
            i += 1
#    print(P1)
        for i5 in range(n):
            P2 += P1[i5]
#            print('P2=', P2)
#        plt.contourf(xi, yi, P2, vmin = 0, vmax = 1, cmap='jet')
        plt.imshow(P2, vmin = P2.min(), vmax = P2.max()/2, cmap='jet', origin='lower', 
           extent=[x.min(), x.max(), y.min(), y.max()])
        plt.pause(0.05)
        if i1==t.min():
            plt.colorbar()
            plt.title('Multipole')
#        print('P1=', P1)
#        print('r=', r)
    return
#        zi1 = scipy.interpolate.griddata((x, y), P1, (xi, yi), method='linear')
    
#    def rho(r, t):
#        rho_plus = (np.sin(0.5*np.pi*(1-r/Ri)))*(np.cos(2*np.pi*f*t))
#        rho_minus = -(np.sin(0.5*np.pi*(1-r/Ri)))*(np.cos(2*np.pi*f*t))
#        return rho_plus, rho_minus
#    fig = plt.figure(1)
#    ax = fig.add_subplot(111)
#    t = 0
#    i = 0
#    r = np.arange(-R_tol, R_tol, itr)
#    return

sound_prop(800, 200, 5, 200, 2, 5)