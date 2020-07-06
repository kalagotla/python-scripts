# Source code for Assignment #2
# Name: YOUR NAME HERE!!!!
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


# Problem #1
def Factorial(x):
    if x < 0 or int(x) != x:
        return('Error: Factorial is defined for non-negative integers only.')
    if x > 0:
        m = 1
        for i in range(1, x+1):
            m = i * m
        return m
    if x == 0:
        return 1
    
# Problem #2
def isNarcissistic(x):
    if x < 0 or int(x) != x:
        return('Error: Narcissistic is defined for non-negative integers only.')
    p = 0
    NumofDigits = len(str(x))
    for i in str(x):
        n = int(i)**NumofDigits
        p = p + n
    if p == x:
        return(True)
    else:
        return(False)
        
# Problem #3
def FindNarcissistic(n):
    x = []
    for m in range(1, 10**9):
        if len(x) <= n:
            p = 0
            NumofDigits = len(str(m))
            for i in str(m):
                i = int(i)**NumofDigits
                p = p + i
            if p == m:
                x += list([m])
            if len(x) == n:
                return x
    
# Problem #4
# Compute the value of pi using the iteration x = x + sin(x)
def ComputePi():
    import math
    x = 3
    y = 3
    for i in range(100):
        x = x + math.sin(x)
        pi = float("{:0.12f}".format(x))
        if y != pi:
            y = pi
        else:
            return pi