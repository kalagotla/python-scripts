# Source code for Assignment #4
# Name: Aatresh Karnam
#
# Required function definitions:
#
#    Problem #1: 
#        Function: Ackley(a, b, c, y)
#        Input: Three scalars (a,b,c) and an array of N elements (y)
#        Output: Tuple (x, f(x))
#
#    Problem #2: 
#        Function: AugLM(N)
#        Input: N (integer)
#        Output: 1D numpy.array containing the solution x
#
def Ackley(a, b, c, y):
    import numpy as np
    import scipy.optimize as opt
    N = np.size(y)    
    x0 = np.zeros(N)
    f = lambda x: a * (1 - np.e**(-b * np.sqrt(sum((x-y)**2)/N))) + np.e - np.e**(sum(np.cos(c*(x-y))/N))
    res = opt.minimize(f, x0, method ="Powell", tol = 1e-6)
    return res,f(res)
    
   
def AugLM(f, g, x0, alpha0, tol):
    import numpy as np
    import scipy.optimize as opt
    alpha = alpha0
    lmbda = np.size(g)
    func = lambda x: f(x) + alpha*g(x)**2 + lmbda*g(x)
    N = np.size(x0)
    delta = np.ones(N)
    while np.linalg.norm(delta. np.inf) > tol:
        residual = opt.minimize(func, x0, method="Powell", tol = 1e-6)
        alpha*= 5
        delta = residual.x - x0
        x0 = residual.x
        lmbda = lmbda + alpha*g(residual.x)
        
        
    return residual, f(residual), g(residual)
    