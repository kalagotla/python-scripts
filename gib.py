# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 20:02:35 2016

@author: mechd
""
# -*- coding: utf-8 -*-
"""
#Spyder Editor

#This is a temporary script file.
import numpy as np
import scipy.optimize as opt
def ForwardEuler(f, g, x0, y0, u0, v0, tRange, N):
    (t, h) = np.linspace(tRange[0], tRange[1], N, retstep=True)
    y = np.zeros(N)
    x = np.zeros(N)
    u = np.zeros(N)
    v = np.zeros(N)
    y[0] = y0
    x[0] = x0
    u[0] = u0
    v[0] = v0
    for k in range(N-1):
        u[k+1] = u[k] + h*f(t[k], u[k], v[k])
        v[k+1] = v[k] + h*g(t[k], u[k], v[k])
        y[k+1] = y[k] + h*v[k]
        x[k+1] = x[k] + h*u[k]
    return (t, y, x)
methods={'Forward Euler': ForwardEuler}

def h(x):
    return(1.6*np.exp(-10*((x-1)**2)))
   
def IVP(f, g,  x0, y0, u0, v0, tRange, N, method="Forward Euler"):
    ivp_method = methods[method]    
    return ivp_method(f, g, x0, y0, u0, v0, tRange, N)  
    
def function(V0):
    import math
    import matplotlib.pyplot as plt
    import scipy.interpolate as interp
    import scipy.optimize as opt
    theta = np.pi/2.1
    x0 = 0
    y0 = 0
    v0=V0*math.cos(theta)
    u0=V0*math.sin(theta)
    k=0
    tRange = [0, 5]
    N = 10000
    probList =[[lambda t,u,v: -k*u*((u*u)+(v*v))**(1/2)],
               [lambda t,u,v: -k*v*(((u*u)+(v*v))**(1/2))-9.8]]
    f  = probList[0][0]; g= probList[1] [0];
    for method in methods:        
        (t,y,x) = IVP(f, g, x0, y0, u0, v0, tRange, N, method=method)
        
    a = y - h(x)
    indx = np.where(a<0)[0][1]
    teqn = interp.interp1d(t, a, kind="cubic")
    tstar = opt.bisect(teqn, t[indx-1], t[indx])
    xeqn = interp.interp1d(x, a, kind="cubic")  
    xstar = opt.bisect(xeqn, x[indx-1], x[indx])
    yeqn = interp.interp1d(y, a, kind="cubic")
    ystar = opt.bisect(yeqn, y[indx-1], y[indx])
    #return xstar-3    
    
#Velocity = opt.bisect(function, 1, 200)
#print(Velocity)
    
    
    #plt.axis([0, 5, 0, 5])
    plt.plot(x, y)
    plt.plot(x, h(x))
#    plt.plot(t, a,'b',t,t-t,'k')


