# -*- coding: utf-8 -*-
"""
Created on Mon Feb 29 16:27:51 2016

@author: mechd
"""

def Ackley(a, b, c, y):
    import numpy as np
    import scipy.optimize as opt
    # import numpy.random as rnd
    # y = rnd.rand(10)
    # x = np.array((1, 2, 3, 4))
    # y = np.array((1, 2))
    N = np.size(y)
    x0 = np.zeros(N)
    # A = lambda x: a * (1 - np.e**(-b * np.sqrt((sum(x[1:]-y[1:])**2)/N)))
    # B = lambda x: np.e - np.e**(sum(np.cos(c*(x[1:]-y[1:])))/N)
    # A = (1-np.e**((-b) * np.sqrt((sum((x-y)**2))/N)))
    # B = np.e - np.e**(sum(np.cos(c*(x-y)))/N)
    # f = lambda x: a * (1 - np.e**(-b * np.sqrt((sum(x-y)**2)/N))) + (np.e - np.e**(sum(np.cos(c*(x-y)))/N))
    f = lambda x: a * (1 - np.e**(-b * np.sqrt(sum((x-y)**2)/N))) + np.e - np.e**(sum(np.cos(c*(x-y))/N))
    res = opt.minimize(f, x0, method="Powell", tol=1e-6)
    # res = opt.minimize_scalar(f, method="brent")
    return res.x, f(res.x)
    # return A, B
    
    # f = a * (1 - np.e(-b * sqrt(sum((x-y)**2)/N))) + np.e - np.e**(np.cos(c*(x-y))/N)
