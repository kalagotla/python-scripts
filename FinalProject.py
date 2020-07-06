# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 13:52:53 2016

@author: mechd
"""

def FinalProject(a, b, lamda, psy0, K):
    import numpy as np
    import scipy.integrate as intg
    import scipy.optimize as opt
    import matplotlib.pyplot as plt
    def integrand(x, y):
        return (y)*K(x, y)
        #return psy(y)*K(x,y)
    def psy(x):
        return psy0(x) + lamda*intg.quad(integrand, a, b, args=(x))[0]
    vec_psy = np.vectorize(psy)
    #x = np.array((1,2))
    #return vec_psy(x)
    x = np.linspace(-1, 1, 50)
    plt.plot(x, vec_psy(x))
    plt.show()