# -*- coding: utf-8 -*-
"""
Created on Sat Oct 13 10:34:20 2018

@author: mechd
"""

def bisection(p):
    import numpy as np
    import matplotlib.pyplot as plt
    #constants
    a = 0
    b = 2-(1/np.sqrt(5))
    Tol = 1e-4
    Nmax = 100
    c = (a+b)/2
    n = 0
    #function definition
    if p == 1:
        f = lambda x: (x-np.sqrt(2))**(-1)
    elif p == 2:
        f = lambda x: (x-1)**7
    elif p == 3:
        def f(x):
            if x>=1:
                a = (x-1)**(1/7)
            elif x<1:
                a = -(1-x)**(1/7)
            return a
    
    #Bisection algorithm
    plt.figure()
    while n<=Nmax and abs(f(c)):
        n = n+1
        if f(a)*f(c)<0:
            b = c
        else:
            a = c
        c = (a+b)/2
        plt.plot(n, a, 'ro')
        plt.plot(n, b, 'bo')
#        plt.plot(n, c, 'co')
    print('sol=', c)
    
    return