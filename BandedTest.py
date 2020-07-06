# Compares the speed of regular LU versus banded LU for a pentadiagoanl 
#    FDM matrix

import numpy as np
import scipy.linalg as sla

M = 5       # Number of times to solve each system
N = 10000    # Size of the system

print("Solving a {0}x{0} linear system using regular LU and banded LU {1} times each".format(N, M))
# Create full matrix
A = np.zeros((N,N))
for k in range(N):  # Diagonal entries
    A[k,k] = 6

for k in range(N-1):  # 1st super/sub digonals
    A[k+1, k] = -4
    A[k, k+1] = -4

for k in range(N-2):  # 2nd super/sub diagonals    
    A[k+2, k  ] = 1
    A[k  , k+2] = 1

fullMem = (A.size * 8)/1024/1024   # Size in megabytes
print("\nFull Matrix Storage Cost: {:.3f}MB".format(fullMem))

# Create banded storage matrix
B = np.zeros((5, N))
B[0,2:N  ] =  1
B[1,1:N  ] = -4
B[2, :   ] =  6
B[3,0:N-1] = -4
B[4,0:N-2] =  1

bandMem = (B.size * 8)/1024/1024   # Size in megabytes
print("Banded Matrix Storage Cost: {:.3f}MB ({:.3%} of full matrix)\n".format(bandMem, bandMem/fullMem))

    
# Create RHS vector
b = np.zeros(N)
b[  0] =  3
b[N-1] =  3
b[  1] = -1
b[N-2] = -1

# Exact solution
xe = np.ones(N)

# Solve the full matrix system M times using LU factorization
for j in range(M):
    x = sla.solve(A, b)
    
r = np.dot(A,x) - b
print("Full Matrix LU Solve:")
print("   Residual = {:.3e}".format(sla.norm(r, np.inf)))
print("   Error = {:.3e}\n".format(sla.norm(x-xe), np.inf))

# Solve the full matrix system M times using LU factorization
for j in range(M):
    x = sla.solve_banded((2,2), B, b)
    
r = np.dot(A,x) - b    # Need to use the full matrix here...
print("Banded Matrix LU Solve:")
print("   Residual = {:.3e}".format(sla.norm(r, np.inf)))
print("   Error = {:.3e}".format(sla.norm(x-xe, np.inf)))

