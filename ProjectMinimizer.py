# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 20:16:27 2016

@author: mechd
"""
# This is minimization file
# H is function
# V0 is 1x2 array (v0, theta)
# distance is the function to be minimized
# input: ProjectTest(lambda x: 1.6*np.e**(-10*(x-1)**2), 0, 9.8, 3)
def ProjectMinimizer(H, k, g, target):
    import numpy as np
    import scipy.optimize as opt
    import importlib
    def theta1():
        if k == 0 and target == 3:
            return np.pi/2.1
        if k == 1:
            return np.pi/3
        if k == 0 and target == 1.25:
            return np.pi/2.0999
    def distance(V0):
        def ForwardEuler2(f, g, x0, y0, x10, y10, tRange, N):
            (t, h) = np.linspace(tRange[0], tRange[1], N, retstep=True)
            y = np.zeros(N)
            x = np.zeros(N)
            y1 = np.zeros(N)
            x1 = np.zeros(N)
            y[0] = y0
            x[0] = x0
            y1[0] = y10
            x1[0] = x10
            for i in range(N-1):
                x1[i+1] = x1[i] + h*f(t[i], x1[i], y1[i])
                y1[i+1] = y1[i] + h*g(t[i], x1[i], y1[i])
                y[i+1] = y[i] + h*y1[i]
                x[i+1] = x[i] + h*x1[i]
            return(t, x, y)
        uPrime = lambda t, u, v: -k*u*(u*u + v*v)**0.5
        vPrime = lambda t, u, v: -k*v*(u*u + v*v)**0.5 - g
        x0 = 0
        y0 = 0
        theta = theta1()
        u0 = V0*np.cos(theta)
        v0 = V0*np.sin(theta)
        tRange = [0, 5]
        N = 10000
        def Function(H):
            H1 = np.zeros(N)
            H1[0] = H(x[0])
            for i in range(N-1):
                H1[i+1] = H(x[i+1])
            return H1
        (t, x, y) = ForwardEuler2(uPrime, vPrime, x0, y0, u0, v0, tRange, N)
        H1 = Function(H)
        def impact():
            for i in range(N-1):
                if y[i+1] - H1[i+1] <= 1e-6:
                    return i+1
        impact = impact()
        return target-x[impact]
    if k == 0:
        a = 0; b = 20
    if k == 1:
        a = 200; b = 300
    theta = theta1()
    alpha = opt.brentq(distance, a, b)
    mod = importlib.import_module("ProjectPlot")
    ProjectPlot = getattr(mod, "ProjectPlot")
    ProjectPlot(H, theta, alpha, k, g, target)
    return alpha, theta, distance(alpha)
    
    
    #return res.x
    #return impact
    #return(t, x, y, H1)
    #plt.plot(x[0:impact], y[0:impact], 'r', x, H1, 'b')
    #plt.plot(x, y, 'r', x, H1, 'b')
    #plt.show()
    
    
    #cons = ({'type': 'ineq', 'fun': lambda V0: V0})
             #{'type': 'ineq', 'fun': lambda V0: np.pi/2-V0[1]},
             #{'type': 'ineq', 'fun': lambda V0: V0[1]},
             #{'type': 'ineq', 'fun': lambda V0: 1e-6-distance(V0)})
    #res = opt.minimize(distance, 10.5, method="cobyla", constraints=cons, options={'tol': 1e-6})
    #res = opt.minimize(distance, (5,np.pi/4), method="L-BFGS-B", bounds=((0, None), (0, np.pi/2)), tol=1e-6)
    #cobyla is better based on inputs. Maybe not considering constraints
    #SLSQP returning same input
    #Nelder-Mead guess is problem
    #res = opt.brute(distance, (0, 100))
    #x0 = 5 # theta = np.pi/2.40101
    #res = opt.minimize(distance, x0, method="Powell", tol=1e-8)
        #100, 400 for k=1, theta = np.pi/2.1, 2.2, 2.3
        # 0, 20 for  k=0, theta = np.pi/3
        # 0, 20 for k=0, theta = np.pi/2.1