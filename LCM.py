# Defines least common multiple (LCM) and greatest common divisor (GCD) functions

def GCD(a, b):
    # Computes the greatest common divisor of two integers
    if (type(a) != int) or (type(b) != int) or (a < 0) or (b < 0):
        print("Error: GCD requires two non-negative integers!")
        return None 
    
    if (a < b):        
        a, b = b, a    # Swap a and b
    
    
    # Euclidean algorithm
    while (b != 0):
        # Set a = b and b = a % b
        a, b = b, (a % b) 
        
    return a


def LCM(a, b):
    # Computes the least common multiple of two integers

    # Input validation
    if (type(a) != int) or (type(b) != int) or (a < 0) or (b < 0):
        print("Error: LCM requires two non-negative integers!")
        return None

    if (a*b == 0):
        return 0 # One or both inputs are zero

    return (a//GCD(a,b))*b


def TestGCD(test_lst=[(2,4,2),(3,11,1),(11,17,1),(6,40,2),(462,1071,21), \
                      (1008,840,168)]):
    # Tests the GCD function for a list of tuples (a, b, GCD(a,b))
    flag = True
    for pair in test_lst:
        a = pair[0]
        b = pair[1]

        gcd = pair[2]
        cGCD = GCD(a,b)

        if (cGCD != gcd):
            print("Error:", cGCD, "!=", gcd, "for a =", a, "and b =", b)
            flag = False
            break

    if (flag == True):
        print("GCD passed tests!")
    else:
       print("GCD failed tests!")
       
def TestLCM(test_lst=[(2,4,4),(3,11,33),(6,40,120),(11,17,187), \
                      (462,1071,23562),(1008,840,5040), (0,2,0), \
                      (2,0,0), (0, 0, 0)]):
    # Tests the LCM function for a list of tuples (a, b, LCM(a,b))
    flag = True
    for pair in test_lst:
        a = pair[0]
        b = pair[1]

        lcm = pair[2]
        cLCM = LCM(a,b)

        if (cLCM != lcm):
            print("Error:", cLCM, "!=", lcm, "for a =", a, "and b =", b)
            flag = False
            break

    if (flag == True):
        print("LCM passed tests!")
    else:
       print("LCM failed tests!")
       

