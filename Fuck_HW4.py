# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 01:35:16 2016

@author: mechd
"""

def Ackley(a, b, c, y):
    import numpy as np
    import scipy.optimize as opt
    N = np.size(y)
    x0 = np.zeros(N)
    f = lambda x: a * (1 - np.e**(-b * np.sqrt(sum((x-y)**2)/N))) + np.e - np.e**(sum(np.cos(c*(x-y))/N))
    res = opt.minimize(f, x0, method="Powell", tol=1e-6)
    return res.x, f(res.x)

def AugLM(f, g, x0, alpha0, tol):
    import numpy as np
    import scipy.optimize as opt
    alpha = alpha0
    fun = lambda x: f(x) + alpha*g(x)**2 + lambda1*g(x)
    N = np.size(x0)
    M = np.size(g)
    delta = np.ones(N)
    lambda1 = np.ones(M)
    while np.linalg.norm(delta, np.inf) > tol:
        res = opt.minimize(fun, x0, method="Nelder-Mead", tol=tol)
        alpha *= 7
        delta = res.x - x0
        lambda1 = lambda1 + alpha*g(res.x)
        x0 = res.x
    return res.x, f(res.x), g(res.x)
