# Source code for Assignment #3
# Name: YOUR NAME HERE!!!!
#
# Required function definitions:
#
#    Problem #1: 
#        Function: LinearLS(x,y)
#        Input: Two 1D numpy.arrays of length N
#        Output: 1D numpy.array of length two [b, m])
#
#    Problem #2: 
#        Function: Problem2(N)
#        Input: Size of the linear system N (integer)
#        Output: 1D numpy.array of length N
#
def LinearLS(x, y):
    import numpy as np
    import numpy.linalg as nla
    # x input vector of size 1xN
    # y input vector of size 1xN
    N = np.size(x)
    M = np.zeros((N, 2))
    for i in range(N):
        M[i, 0] = 1
        M[i, 1] = x[i]
    A = np.dot(np.transpose(M), M)
    B = np.dot(np.transpose(M), np.transpose(y))
    # solving Ax = B
    for i in range(5):
        z = nla.solve(A, B)
    return z


def SolveFDM(N):
    import numpy as np
    import scipy.linalg as sla
    # Matrices A and b
    A = np.zeros((N, N))
    b = np.zeros(N)
    for i in range(N):
        b[i] = - 2 / (N + 1)**2
        for j in range(N):
            if i == j:
                A[i, j] = 2
            elif i+1 == j or i-1 == j:
                A[i, j] = -1
                A[i, j] = -1
    B = np.zeros((3, N))
    B[0, 1:N] = -1
    B[1, :] = 2
    B[2, 0:N-1] = -1
    # Banded matrix solver
    for i in range(5):
        x = sla.solve_banded((1, 1), B,b)
    return x
