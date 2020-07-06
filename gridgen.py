# -*- coding: utf-8 -*-
"""
Created on Sat Jun 10 12:15:18 2017

@author: mechd
"""

def grid(theta, d_theta, n, tol):
    """"theta = pi/2 - d_theta*n*tol"""
#    grid(0.942478, 2.51327, 25, 0.01)
    import numpy as np
    import math
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    import random
    a = np.zeros((n*n*n, 3))
    b = np.zeros((n*n*n, 3))
    c = np.zeros((2*n*n*n, 3))
    kc = 0
    kd = 0
    kt1 = 0
    kt2 = 0
    for k in range(n):
        for j in range(n):
            m = 0
            for i in range(n):
                m = m + tol
#            print(m)
                a[kc, 0] = (i-1) + (j-1)/np.tan((np.pi/2)-m*d_theta)
                a[kc, 1] = j
                a[kc, 2] = k
                kc = kc+1
    for k in range(n):
        for j in range(n):
            m = 0
            for i in range(n):
                m = m + tol
#            print(m)
                b[kd, 0] = (i+n-1) + (j-1)/np.tan(theta+m*d_theta)
                b[kd, 1] = j
                b[kd, 2] = k
                kd = kd+1
    for k in range(n):
        for j in range(n):
            for i in range(n):
                c[kt1, 0] = a[kt2, 0]
                c[kt1, 1] = a[kt2, 1]
                c[kt1, 2] = a[kt2, 2]
                c[kt1+n, 0] = b[kt2, 0]
                c[kt1+n, 1] = b[kt2, 1]
                c[kt1+n, 2] = b[kt2, 2]
                kt2 = kt2+1
                kt1 = kt1+1
            kt1 = kt1+n
#    print(a)
#    plt.plot(a,b,'ro')
#    plt.plot(c,d,'bo')
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(a[0:n**3-1, 0], a[0:n**3-1, 1], a[0:n**3-1, 2], c='r')
    ax.scatter(b[0:n**3-1, 0], b[0:n**3-1, 1], b[0:n**3-1, 2], c='b')
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    np.savetxt('griddata.o', c/10000)
    
    
    return