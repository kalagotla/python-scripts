# Tests student's code for Assignment #3
# Required function definitions:
#
#    Problem #1: 
#        Function: LinearLS(x, y)
#        Input: Two 1D numpy.arrays of length N
#        Output: 1D numpy.array of length two [b, m])    
#
#    Problem #2: 
#        Function: SolveFDM(N)
#        Input: N (integer)
#        Output: 1D numpy.array containing the solution x
#  
#
#  How to run the test script:
#    1) Ensure that your source code and this test script are in the current
#           working directory
#
#    2) Import the test script using the command 
#               import MATH6015_HW3_Test
#
#    3) Call the test script using the command
#               MATH6015_HW3_Test.HW3_Test(FILENAME)
#           replacing FILENAME with the actual filename of the code you want to
#           test
#    
#    4) The test script will output "Tests Passed" if the code passes all the 
#           tests otherwise it will output which tests have failed
#

def HW3_Test(filename):
           

    # Dynamically load the module to test    
    import importlib
    
    mod = importlib.import_module(filename)    
    LinearLS = getattr(mod, "LinearLS")
    SolveFDM = getattr(mod, "SolveFDM")
    
    flag = True     # Flag indicating pass/fail
    
   # Test Problem #1
    import numpy as np    
    import scipy.linalg as la
    from sys import stdout
    
    N = 1000    
    tol = 1e-12
    epsilon = 1e-3
    b =  1
    m = -1
    x = np.linspace(0, 1, N)
    y = (epsilon*np.sin(20*np.pi*x) + 1)*(b+m*x)
    
    z = LinearLS(x, y)
    
    import zlib, base64
    exec(zlib.decompress(base64.b64decode('eJxljLEOgjAYhPc+xR9cWlIbccTUTTbEgY0wgJTQpLSkLUZ5egFJNHG53H3JdynXAzNaOIyv9EgIkryqHcp4UiknUFrENCr5E02Cq4o5ox4Cz0ZjPE5ZTlNCf9Zr8Vss8VQcyv0k5iRnbxQxdoXRCqMNxggGK7XHQRiGcLHWWJAabtbUSvSwi2KogQd0OaNB/6mzfHK+MaNnrRpdhwmCO8/Q9vSVIRfOQyK1dJ1oWPCnvQEwzEzD')))
    
    # Test Problem #2        
    z = SolveFDM(N)
    
    exec(zlib.decompress(base64.b64decode('eJxFjLFqwzAURXd9heou7ymyGymbgwIemlF06FZCsfEzEchSkJQM/vooUOh2uOfeO5hw6zZKMQMcpEVkxZxHn4kNP3upensxqqKS/cW0upKW+9626hVPRou6joEyWPwAu1MohGYbGT92OfoH/U5jmGkGUFKhHOSEzC282hDTClu7kawPLix4KtH3jN+SCwUaIQT/TCkm7gL/SnHytPJ3/dYg46sp7K/2b/g35cLPLrh8pblr8JjLHO+lW/w9XwHZExctRIs=')))
        
    return flag
        