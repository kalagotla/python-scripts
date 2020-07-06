# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 14:20:00 2016

@author: mechd
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 20:16:27 2016

@author: mechd
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 14:09:08 2016

@author: mechd
"""
# This is minimization file
# H is function
# V0 is 1x2 array (v0, theta)
# distance is the function to be minimized
# input: ProjectTest(lambda x: 1.6*np.e**(-10*(x-1)**2), 0, 9.8, 3)
def ProjectTest(H, k, g, target):
    import numpy as np
    import scipy.integrate as intg
    import scipy.optimize as opt
    import matplotlib.pyplot as plt
    import importlib
    mod = importlib.import_module("ProjectPlot")
    ProjectPlot = getattr(mod, "ProjectPlot")
    def theta1():
        if k == 0 and target == 3:
            return np.pi/2.1
        if k == 1:
            return np.pi/3
        if k == 0 and target == 1.25:
            return np.pi/2.0999
    def distance(V0):
        theta = theta1()
        ProjectPlot(H, theta, alpha, k, g, target)
        return target-x[imapct]
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
    if k == 0:
        a = 0; b = 20
    if k == 1:
        a = 200; b = 400
    theta = theta1()
    alpha = opt.brentq(distance, a, b)
    #ProjectPlot(H, theta, alpha, k, g, target)
    return alpha, theta, distance(alpha)
    #return res.x
    #return impact
    #return(t, x, y, H1)
    #plt.plot(x[0:impact], y[0:impact], 'r', x, H1, 'b')
    #plt.plot(x, y, 'r', x, H1, 'b')
    #plt.show()