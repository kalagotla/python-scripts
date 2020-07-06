# Tests student's code for Assignment #2
# Required function definitions:
#
#    Problem #1: 
#        Function: Factorial(n)
#        Input: n (integer)
#        Output: n! (integer)       
#
#    Problem #2: 
#        Function: isNarcissistic(n)
#        Input: n (integer)
#        Output: True/False (boolean)
#  
#    Problem #3:
#        Function: FindNarcissistic(N)
#        Input: N (integer) # of narcissistic numbers to find
#        Output: List of the first N narcissistic numbers
#
#    Problem #4:
#        Function: ComputePi()
#        Input: None
#        Ouput: Value of pi to at least 12 digits of accuracy
#
#
#  How to run the test script:
#    1) Ensure that your source code and this test script are in the current
#           working directory
#
#    2) Import the test script using the command 
#               import MATH6015_HW2_Test
#
#    3) Call the test script using the command
#               MATH6015_HW2_Test.HW2_Test(FILENAME)
#           replacing FILENAME with the actual filename of the code you want to
#           test
#    
#    4) The test script will output "Tests Passed" if the code passes all the 
#           tests otherwise it will output which tests have failed
#
#    5) Depending on how fast your computer is, Test #3 might take a few
#           seconds to a minute to finish
#   
#    6) If the test script hangs, then you probably have an error in the 
#           FindNarcissistic function and the script cannot find the required 
#           twenty numbers

def HW2_Test(filename):
    # Dynamically load the module to test
    import importlib
    
    mod = importlib.import_module(filename)
    Factorial = getattr(mod, "Factorial")
    isNarcissistic = getattr(mod, "isNarcissistic")
    FindNarcissistic = getattr(mod, "FindNarcissistic")
    ComputePi = getattr(mod, "ComputePi")
    
    from sys import stdout
    
    flag = True     # Flag indicating pass/fail
    
    # Test Problem #1
    import math
    n_lst = list(range(11)) + [20, 30, 40]
    for n in n_lst:
        f1 = math.factorial(n)
        f2 = Factorial(n)
        if (f1 - f2) != 0:
            print("*** Error in Problem #1: n =", n, "Factorial(n) =", f2, "factorial(n) =", f1); stdout.flush()
            flag = False
    print("Problem #1 Test Finished."); stdout.flush()
    
    # Test Problem #2
    true_lst = list(range(1,10)) + [153, 9474, 93084, 32164049651]
    false_lst = [17, 186, 1395, 65479, 984315, 6527416, 23486475]
    for n in true_lst:
        if (not isNarcissistic(n)):
            print("*** Error in Problem #2: n =", n, "is narcissitic but function returned false"); stdout.flush()
            flag = False
    for n in false_lst:
        if (isNarcissistic(n)):
            print("*** Error in Problem #2: n =", n, "is not narcissitic but function returned true"); stdout.flush()
            flag = False
    print("Problem #2 Test Finished."); stdout.flush()        
    
    # Test Problem #3
    narc_lst = [1, 2, 3, 4, 5, 6, 7, 8, 9, 153, 370, 371, 407, 1634, 8208, 9474, 54748, 92727, 93084, 548834]
    test_lst = FindNarcissistic(20)
    if (narc_lst != test_lst):
        print("*** Error in Problem #3: Compute list of narcissistic numbers is not correct")
        print("***    ", narc_lst)
        print("***    ", test_lst); stdout.flush()
        flag = False
    print("Problem #3 Test Finished."); stdout.flush()            
            
    # Test Problem #4
    import math
    cPi = ComputePi()
    delta = abs(math.pi - cPi)
    if (delta > 1e-12):
        print("*** Error in Problem #4: Estimate of pi is not within 1e-12.")
        print("***     Computed pi =", cPi, " delta =", delta); stdout.flush()
        flag = False
    print("Problem #4 Test Finished."); stdout.flush()            
    
    if flag:
        print("Tests passed!")
    else:
        print("Tests failed!")
        
    return flag
        