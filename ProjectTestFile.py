# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 17:35:14 2016

@author: mechd
"""
def ProjectTestFile(filename):
    import importlib
    import numpy as np
    import scipy.integrate as intg
    import scipy.optimize as opt
    import matplotlib.pyplot as plt
    
    mod = importlib.import_module(filename)
    ProjectTest2 = getattr(mod, "ProjectTest2")
    H1 = lambda x: 1.6*np.e**(-10*(x-1)**2)
    k1 = 0
    g1 = 9.8
    target1 = 3
    f = ProjectTest2(H1, k1, g1, target1)
    H2 = lambda x: 1.6*np.e**(-10*(x-1)**2)
    k2 = 1
    g2 = 9.8
    target2 = 3
    g = ProjectTest2(H2, k2, g2, target2)
    H3 = lambda x: 3*np.e**(-10*(x-1)**2)
    k3 = 0
    g3 = 9.8
    target3 = 1.25
    h = ProjectTest2(H3, k3, g3, target3)
    return f, g, h
    #Prblm 1: ProjectTest2(lambda x: 1.6*np.e**(-10*(x-1)**2), 0, 9.8, 3)
    #Prblm 2: ProjectTest2(lambda x: 1.6*np.e**(-10*(x-1)**2), 1, 9.8, 3)
    #Prblm 3: ProjectTest2(lambda x: 3*np.e**(-10*(x-1)**2), 0, 9.8, 1.25)
    