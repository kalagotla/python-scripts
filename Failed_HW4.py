# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 09:12:52 2016

@author: mechd
"""

# Source code for Assignment #4
# Name: Dilip Kalagotla
#
# Required function definitions:
#
#    Problem #1: 
#        Function: Ackley(a, b, c, y)
#        Input: Three scalars (a,b,c) and an array of N elements (y)
#        Output: Tuple (x, f(x))
#
#    Problem #2: 
#        Function: AugLM(N)
#        Input: N (integer)
#        Output: 1D numpy.array containing the solution x
#
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

def AugLM(f, g, x0, alpha0, tol):
    import numpy as np
    import scipy.optimize as opt
    alpha = alpha0
    fun = lambda x: f(x) + alpha*g(x)**2
    N = np.size(x0)
    delta = np.ones(N)
    while np.linalg.norm(delta, np.inf) > tol:
        res = opt.minimize(fun, x0, method="Nelder-Mead", tol=tol)
        alpha *= 7
        delta = res.x - x0
        x0 = res.x
    return res.x, f(res.x), g(res.x)
