# --------------------------------------------------------------------------- #
def SolveMatrixSystem(A, C):
# Solve the matrix system AB = C for B given A and C
# This function can compute the inverse of A if C is the identity matrix
#   but numpy.linalg.inv is faster
#
# A [in]: NxN numpy array
# C [in]: NxN numpy array
# 
# B [out]: NxN numpy array

# Insert input checking here

    import numpy as np
    import scipy.linalg as la
    
    (N, N) = C.shape
    
    (LU, P) = la.lu_factor(A)
    
    B = np.zeros((N,N))
    
    for k in range(N):
        c = C[k,:]
        z = la.lu_solve((LU, P), c)
        B[k,:] = z
        
    return B
# --------------------------------------------------------------------------- #
    
# --------------------------------------------------------------------------- #
def SolveBlockTriDiag(A, B, C, b):
# Solve the block tridiagonal matrix system Mx = b where the diagonal blocks 
#   of M are given by A, the superdiagonal blocks by B, and the subdiagonal 
#   blocks by C
#
# A [in]: RxRxM numpy array (M RxR matrices)
# B [in]: RxRxM-1 numpy array (M-1 RxR matrices)
# C [in]: RxRxM-1 numpy array (M-1 RxR matrices)
# b [in]: 1D numpy array of length N = MR
#
# x [in]: 1D numpy array of length N = MR 
#

    import numpy as np
    import scipy.linalg as la
    
    # Insert input checking here
    
    (N,) = b.shape
    (R, R, M) = A.shape
    
    U = np.zeros((R, R, M))
    L = np.zeros((R, R, M-1))
    
    # Block Tridiagonal LU Factorization
    U[:, :, 0] = A[:, :, 0]
    for k in range(M-1):
        L[:, :, k  ] = np.dot(la.inv(U[:, :, k]), C[:,:,k])
        U[:, :, k+1] = A[:, :, k+1] - np.dot(L[:, :, k], B[:, :, k])
    
    # Forward Solve
    y = np.copy(b)    
    for k in range(M-1):
        idx1 =  k   *R
        idx2 = (k+1)*R
        idx3 = (k+2)*R        
        y[idx2:idx3] = b[idx2:idx3] - np.dot(L[:, :, k], y[idx1:idx2])
        
    # Backsolve
    x = np.zeros(N)
    idx1 = N-R
    idx2 = N
    x[idx1:idx2] = la.solve(U[:, :, M-1], y[idx1:idx2])
    for k in range(M-2, -1, -1):
        idx1 =  k   *R
        idx2 = (k+1)*R
        idx3 = (k+2)*R        
        x[idx1:idx2] = la.solve(U[:, :, k], y[idx1:idx2] - np.dot(B[:, :, k], x[idx2:idx3]))
        
    return x    

# --------------------------------------------------------------------------- #


# --------------------------------------------------------------------------- #
def CreateFullMatrixFromBlocks(A, B, C):
    import numpy as np
    (p, p, m) = A.shape
    N = p*m
    M = np.zeros((N, N))
        
    M[N-m:N, N-m:N] = A[:, :, m-1]
    for k in range(m-1):
        idx1 =  k   *m
        idx2 = (k+1)*m
        idx3 = (k+2)*m
        M[idx1:idx2, idx1:idx2] = A[:, :, k]
        M[idx1:idx2, idx2:idx3] = B[:, :, k]
        M[idx2:idx3, idx1:idx2] = C[:, :, k]
        
    return M
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
def TestTriBlockSolver(m):
    import numpy as np    
    import numpy.linalg as la
        
    N = m**2
    
    A = np.zeros((m, m, m  ))
    B = np.zeros((m, m, m-1))
    #C = np.zeros((m, m, m-1)) # Not needed since B = C = I
    
    Afull = np.zeros((N,N))
    
    for k in range(m):
        A[:, :, k] = np.diagflat(4*np.ones(m), 0)
        
    for k in range(m-1):
        B[:, :, k] = np.eye(m)
        
    Afull = CreateFullMatrixFromBlocks(A, B, B)
    
    xe = np.random.rand(N)
    b = np.dot(Afull, xe)
    
    x = SolveBlockTriDiag(A, B, B, b)
    x2 = la.solve(Afull, b)
    
    delta = la.norm(x-xe, np.inf)
    print("Absolute Error:", str(delta))
    print("Relative Rrror:", str(delta/la.norm(xe,np.inf)))
    
    delta = la.norm(x-x2, np.inf)
    print("Absolute Delta:", str(delta))
    print("Relative Delta:", str(delta/la.norm(x2,np.inf)))    
        
# --------------------------------------------------------------------------- #
    
    