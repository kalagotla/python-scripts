# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 13:28:08 2017

@author: mechd
"""

def flow(vf1, vf2, rho1, rho2, T1, T2, delta, theta, d_theta, n, tol):
#    flow(521, 428, 0.161443, 0.942478, 2.51327, 25, 0.01)
    import numpy as np
    import math
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    import random
#    OSWPI(1.5, 1.273, 300, 1.4, 286.9, 0.161443, 1000000, 1e-3)
#    grid(0.942478, 2.51327, 25, 0.01)
    V1 = np.zeros((n*n*n, 5))
    V2 = np.zeros((n*n*n, 5))
    V2t = np.zeros(n*n*n)
    V = np.zeros((2*n*n*n, 5))
    kc = 0
    kd = 0
    kt1 = 0
    kt2 = 0
    for k in range(n):
        for j in range(n):
            for i in range(n):
                V1[kc, 0] = rho1
                V1[kc, 1] = vf1
                V1[kc, 2] = 0
                V1[kc, 3] = 0
                V1[kc, 4] = T1
                kc = kc+1
    for k in range(n):
        for j in range(n):
            m = 0
            for i in range(n):
                m = m + tol
                V2[kd, 0] = rho2
                V2[kd, 1] = vf2*np.cos(delta)
                V2[kd, 2] = vf2*np.sin(delta)
#                V2[kd, 1] = vf2*(np.cos(delta) - np.sin(delta)/np.tan(theta+m*d_theta))
#                V2[kd, 2] = vf2*np.sin(delta)/np.sin(theta+m*d_theta)
                V2[kd, 3] = 0
                V2[kd, 4] = T2
                V2t[kd] = np.sqrt(V2[kd, 0]**2 + V2[kd, 1]**2 + 2*V2[kd, 0]*V2[kd, 1]*np.cos(theta+m*d_theta))
                kd = kd+1
    for k in range(n):
        for j in range(n):
            for i in range(n):
                V[kt1, 0] = V1[kt2, 0]
                V[kt1, 1] = V1[kt2, 1]
                V[kt1, 2] = V1[kt2, 2]
                V[kt1, 3] = V1[kt2, 3]
                V[kt1, 4] = V1[kt2, 4]
                V[kt1+n, 0] = V2[kt2, 0]
                V[kt1+n, 1] = V2[kt2, 1]
                V[kt1+n, 2] = V2[kt2, 2]
                V[kt1+n, 3] = V2[kt2, 3]
                V[kt1+n, 4] = V2[kt2, 4]
                kt2 = kt2+1
                kt1 = kt1+1
            kt1 = kt1+n
#    print(V[2*15612:2*15624])
    np.savetxt('flowdata.o', V)
    return