# -*- coding: utf-8 -*-
"""
Created on Sun Feb 21 12:35:36 2016

@author: mechd
"""

def LinearLS(x, y):
    import numpy as np
    import numpy.linalg as nla
    # x input vector of size 1xN
    # y input vector of size 1xN
    # x = np.matrix('1, 2, 3')
    # y = np.matrix('2, 4, 6')
    N = np.size(x)
    M = np.zeros((N, 2))
    for i in range(N):
        M[i, 0] = 1
        M[i, 1] = x[0, i]
    A = np.dot(np.transpose(M), M)
    B = np.dot(np.transpose(M), np.transpose(y))
    # solving Ax = B
    for i in range(5):
        z = nla.solve(A, B)
    return z
