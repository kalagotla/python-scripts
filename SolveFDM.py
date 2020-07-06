# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 14:40:59 2016

@author: mechd
"""

def SolveFDM(N):
    import numpy as np
    import scipy.linalg as sla
    import numpy.linalg as nla
    A = np.zeros((N, N))
    b = np.zeros((N, 1))
    for i in range(N):
        b[i, 0] = 2 / (i + 2)**2
        for j in range(N):
            if i == j:
                A[i, j] = 2
            elif i+1 == j or i-1 == j:
                A[i, j] = -1
                A[i, j] = -1
    # This is LU factorization           
    # return A
    # u = np.zeros((N, N))
    # l = np.zeros((N, N))
    # for k in range(1, N):
    #    u[0, 0] = L[0, 0]
    #   l[k, k-1] = L[k, k-1] / u[k-1, k-1]
    #    u[k, k] = L[k, k] - l[k, k-1] * L[k-1, k]
    # return l
    # Banded storage matrix
    B = np.zeros((3, N))
    B[0, 1:N] = -1
    B[1, :] = 2
    B[2, 0:N-1] = -1
    # return B
    x = sla.solve_banded((1, 1), B, b)
    z = nla.solve(A, b)
    return x, z
    