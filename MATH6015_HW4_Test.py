# Tests student's code for Assignment #4
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
#
#  How to run the test script:
#    1) Ensure that your source code and this test script are in the current
#           working directory
#
#    2) Import the test script using the command 
#               import MATH6015_HW4_Test
#
#    3) Call the test script using the command
#               MATH6015_HW4_Test.HW4_Test(FILENAME)
#           replacing FILENAME with the actual filename of the code you want to
#           test
#    
#    4) The test script will output "Tests Passed" if the code passes all the 
#           tests otherwise it will output which tests have failed
#

def HW4_Test(filename):
           

    # Dynamically load the module to test    
    import importlib
    
    mod = importlib.import_module(filename)    
    Ackley = getattr(mod, "Ackley")
    AugLM = getattr(mod, "AugLM")
    
    flag = True     # Flag indicating pass/fail
    
    import numpy as np        
    import numpy.random as rnd
    import scipy.optimize as opt
    from sys import stdout
    import zlib, base64
            
   # Test Problem #1
    tol = 1e-6
    a = 20
    b = 0.2
    c = 2*np.pi
    
    N = 10
    y = rnd.rand(N)
    
    exec(zlib.decompress(base64.b64decode('eJxljj0LwjAURff8imddkpIWXJUILpkdHFxT+2KDrx80qaSQHy9aRcHx3Mu5XB6lvRsS6nC5Ec7cyEpe5CxYjRTMWXVD2ZrITeV5LGbxzrV6Bi+RDaPrAs9SilDAGVMClcnFllmyX9SCOcuXZh96Ev24kH7RlsF7Ks9zOI59RdjCegPaOMIaQoNwQh9Wmdj5UPdTKC1NvuGCgSVzVdqQx8+dH/8pgXad8w3W5b/9AMh9VGY=')))
    
    # Test Problem #2         
    g = []
    g.append(lambda x: x[0] - x[1])
    g.append(lambda x: x[0]**2 + x[1]**2 - 1)
    g.append(lambda x: (x[0] - 0.5)**2 + x[1]**2 - 0.125)
    g.append(lambda x: np.dot(x,x) - 1)
    
    N = 2    
    x0 = np.ones(N)
    alpha0 = 10
    
    for f in g:
        (x, fval, gval) = AugLM(opt.rosen, f, x0, alpha0, tol)              
        exec(zlib.decompress(base64.b64decode('eJw1zDFPhDAYxvG9n6LWoS0pxDhi6niTyQHnYGIcevL2aK60HH0xqPjd5Yhuv//wPHsdByx6F1zvvkBcY4wJgprvVA/YxVazw9Ohrph6jyHhaFzApF+/OX4OwEsOF664nQIv7c+bJLU2x0TOemd8AlLpMBS9mUUt9sWcz1KSRl+9DnL7YbwkuPZp0zCu14Ity0xz+gLLQjVTlWKL3dSsOq2iTKEkzorqEaOXcRTNP3BDSejfVZZltBrj0UNPb+/pzjgPLcUO6DMkvGHyIWEbJyysn1InJKEXfSa/QP1W8w==')))

        
    print("Problem #2 Test Finished.")
    
    return flag
        