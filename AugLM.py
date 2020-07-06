# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 08:25:18 2016

@author: mechd
"""

def AugLM(alpha0, tol):
    import numpy as np
    import scipy.optimize as opt
    f = lambda x: x[0]**2 + x[1]**2 - 1
    g = lambda x: x[0] - x[1]
    fun = lambda x: f(x) + alpha*g(x)**2
    guess = np.zeros(2)
    delta = np.ones(2)
    alpha = alpha0
    while np.linalg.norm(delta, np.inf) > tol:
        res = opt.minimize(fun, guess, method="Nelder-Mead", tol=tol)
        alpha *= 4
        delta = res.x - guess
        guess = res.x
    return res.x, f(res.x)