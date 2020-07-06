# Source code for Assignment #4
# Name: YOUR NAME HERE!!!!
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
    # y = np.array((1, 2, 3, 4))
    # y = np.array((1, 2))
    N = np.size(y)
    guess = np.zeros(N)
    delta = np.ones(N)
    # A = lambda x: a * (1 - np.e**(-b * np.sqrt((sum(x[1:]-y[1:])**2)/N)))
    # B = lambda x: np.e - np.e**(sum(np.cos(c*(x[1:]-y[1:])))/N)
    # A = (a) * (1-np.e**((-b) * np.sqrt((sum((x-y)**2))/N)))
    # B = (np.e**(1) - np.e**(sum(np.cos((c) * (x-y))))/N)
    f = lambda x: a * (1 - np.e**(-b * np.sqrt((sum(x-y)**2)/N))) + np.e - np.e**(sum(np.cos(c*(x-y)))/N)
    while np.linalg.norm(delta, np.inf) > 1e-6:
        res = opt.minimize(f, guess, method="Nelder-Mead", tol=1e-6)
        delta = res.x - guess
        guess = res.x
    # res = opt.minimize_scalar(f, method="brent")
    #print(res.x, sum(res.x))
    return res.x, f(res.x)
    
   
def AugLM(f, g, x0, alpha0, tol):
    import numpy as np
    import scipy.optimize as opt
    fun = lambda x: f(x) + alpha*g(x)**2
    N = np.size(x0)
    delta = np.ones(N)
    alpha = alpha0
    while np.linalg.norm(delta, np.inf) > tol:
        res = opt.minimize(fun, x0, method="Nelder-Mead", tol=tol)
        alpha *= 4
        delta = res.x - x0
        x0 = res.x
    return res.x, f(res.x), g(res.x)
