# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 14:09:08 2016

@author: mechd
"""
# input: Project(lambda x: 1.6*np.e**(-10*(x-1)**2), np.pi/3, 10, 0, 9.8, 3)
def Project(H, theta, V0, k, g, target):
    import numpy as np
    import scipy.integrate as intg
    import scipy.optimize as opt
    import matplotlib.pyplot as plt
    def ForwardEuler2(F, G, x0, y0, x10, y10, tRange, N):
        (t, h) = np.linspace(tRange[0], tRange[1], N, retstep=True)
        y = np.zeros(N)
        x = np.zeros(N)
        y1 = np.zeros(N)
        x1 = np.zeros(N)
        y[0] = y0
        x[0] = x0
        y1[0] = y10
        x1[0] = x10
        for i in range(N-1):
            x1[i+1] = x1[i] + h*F(t[i], x1[i], y1[i])
            y1[i+1] = y1[i] + h*G(t[i], x1[i], y1[i])
            y[i+1] = y[i] + h*y1[i]
            x[i+1] = x[i] + h*x1[i]
        return (t, x, y, x1, y1)
    uPrime = lambda t, u, v: -k*u*((u*u + v*v)**0.5)
    vPrime = lambda t, u, v: -k*v*((u*u + v*v)**0.5) - g
    x0 = 0
    y0 = 0
    u0 = V0*np.cos(theta)
    v0 = V0*np.sin(theta)
    tRange = [0, 10]
    N = 100000
    def Function(H):
        H1 = np.zeros(N)
        H1[0] = H(x[0])
        for i in range(N-1):
            H1[i+1] = H(x[i+1])
        return H1
    (t, x, y, x1, y1) = ForwardEuler2(uPrime, vPrime, x0, y0, u0, v0, tRange, N)
    #return x[0:1000]
    H1 = Function(H)
    #return np.size(H1), np.size(x), x
    def impact():
        for i in range(N-1):
            if y[i+1] - H1[i+1] <= 1e-12:
                return i+1
    impact = impact()
    #return impact
    distance = x[impact]-target
    #return(t, x, y, H1)
    #return distance
    plt.figure()
    plt.plot(x[0:impact], y[0:impact], 'r', x, H1, 'b', target, 0, '*')
    return x[impact]
    #plt.plot(x[0:100], y[0:100], 'r', x, H1, 'b')
    #plt.plot(x, y, 'r', x, H1, 'b')
    plt.show()