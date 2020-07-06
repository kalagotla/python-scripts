# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 10:04:06 2017

@author: mechd
"""

def Amit(iter):
    import numpy as np
    phi_0 = 0
    phi_n = phi_0
    for i in range(iter):
        dx = 1/iter
        phi_n = phi_n + dx
        print(phi_n)
    return