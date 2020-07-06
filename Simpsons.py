# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 18:25:29 2016

@author: mechd
"""
def FinalProject1(a, b, lamda, psy0, K, problem):
    import numpy as np
    from scipy.interpolate import interp1d
    import matplotlib.pyplot as plt
    import numpy.linalg as la
    psy_max1 = 0
    Error1 = 1
    psy_avg1 = 0
    print('Problem   ', 'N       ', ' psy_avg           ', ' Error         ', '       Convergence')
    for i in (2**x for x in range(6)):
        N = i*100
        h = (b-a)/N
        g = (b-a)/N
        x = np.zeros(N+1)
        y = np.zeros(N+1)
        K1 = np.zeros((N+1, N+1))
        for i in range(N+1):
            for j in range(N+1):
                K1[i, j] = K(a+h*i, a+h*j)
        I = np.eye(N+1, N+1)
        W = np.diag(np.ones(N+1), k=0)
        for i in range(N+1):
            if i != 0 and i != N:
                if i%2 == 0:                
                    W[i, i] = 2
                else:
                    W[i,i] = 4
        W = ((b-a)/(3*N))*W
        P = la.inv(I-lamda*np.dot(K1, W))
        psy_0 = np.zeros(N+1)
        for i in range(N+1):
            psy_0[i] = psy0(a+h*i)
        psy = np.dot(P, psy_0)
        plt.plot(np.linspace(a, b, N+1), psy)
        plt.show()
        psy_inter = interp1d(np.linspace(a, b, N+1), psy, kind='cubic')
        psy1 = psy_inter(a)
        psy2 = psy_inter(a+3*(b-a)/4)
        psy3 = psy_inter(a+(b-a)/4)
        psy4 = psy_inter(a+(b-a)/2)
        psy5 = psy_inter(b)
        psy_avg = (psy1+psy2+psy3+psy4+psy5)/3
        psy_max = la.norm(psy, np.inf)
        #Error = psy_max1 - psy_max
        Error = psy_avg1-psy_avg
        psy_max1 = psy_max
        psy_avg1 = psy_avg
        Convergence = Error1/Error
        Error1 = Error
        if N != 100 and N != 200:
            print(problem, "       ", i, "    ", psy_avg, "    ", Error, "    ", Convergence)
    plt.figure()
FinalProject1(-1, 1, -1, lambda x: np.e**x, lambda x, y: x*np.e**(y*(1-x)), 1)
FinalProject1(0, np.pi, -1, lambda x: np.sin(10*x), lambda x, y: np.sin(x+y), 2)
FinalProject1(0, np.pi, -1, lambda x: 1+np.sin(np.pi*x), lambda x, y: x*np.cos(x*y), 3)