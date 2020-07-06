# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 20:54:10 2016

@author: mechd
"""

#import numpy as np
#import scipy.integrate as intrg
#import matplotlib.pyplot as plt
#t = 100
#u = lambda u, v: -np.sqrt(u**2+v**2)*u
#v = lambda u, v: -np.sqrt(u**2+v**2)*v-9.8
#u0= (np.sqrt(3)/2, .5)
#s =intrg.odeint(u, u0, t)
#print(s)
#
#y = lambda x: 3*np.exp(-10*(x-1)**2)
#plt.plot(t, y)
#plt.show()


#import numpy as np
#import scipy.integrate as intrg
#import matplotlib.pyplot as plot
#import math as math
#
#N=100
#b = 10
#
## timestep on x axis
#t = np.linspace(0, 5, 100)
#
## height function
#y = np.exp(-10*(t-1)**2)*1.6
#
##projectile motion x' and y' (solved 2nd order by hand)
#xp = lambda x,y: -np.sqrt(x**2+y**2)*x
#yp = lambda y,x: -np.sqrt(x**2+y**2)*y - 9.8
#
#
#initialconditions_xp = (45*np.sqrt(2), 45*np.sqrt(2))
#initialconditions_yp = (90*np.sqrt(2), 90*np.sqrt(2))
#
## ode solver for x and y
#s = intrg.odeint(xp , initialconditions_xp, t)
#u = intrg.odeint(yp, initialconditions_yp, t)
#
## plotting 
#plot.plot(t, y, 'b', t, s, 'r', t, u, 'g')
#
#plot.show()
import numpy as np
import scipy.integrate as intr
import scipy.interpolate as interp
import matplotlib.pyplot as plt
import scipy.optimize as opt

def f(y, t):  #defining ODEs
    u = y[2]
    v = y[3]
    g = 9.8
    k=0  # or 1
    fv = np.zeros(4)
    fv[0] = u
    fv[1] = v
    fv[2] = -k*np.sqrt(u**2+v**2)*u
    fv[3] = -k*np.sqrt(u**2+v**2)*v-g
    return fv
v0 = 5
theta = np.pi/3   
y0= np.array([0, 0, v0*np.cos(theta), v0*np.sin(theta)])  #initial conditions
t = np.linspace(0, 1, 1000) # timestep
sol = intr.odeint(f, y0, t)
x = sol[:,0]
y = sol[:,1]
def H(x):  #the gound function
    return 1.6*np.exp(-10*(x-1)**2)
#print(H(x))
a = y - H(x)   # 
#print(np.size(a))
A = interp.interp1d(t, a, kind='linear')
idx = np.where(a<=0)[0]
idy = np.where(a>=0)[0]
#print(idx)
intX = a[idx][1]
intZ = a[idy][1]
print(intX, intZ)
#A = interp.interp1d(t, a, kind='cubic')
#intX = interp.interpolate(x, y, kind='cubic')
#
r = opt.bisect(A, intX, intZ)
#print(r) 
#
#z = np.array([theta], [v0])
#def compd2(z):
#    
#    return d2
    
#plt.plot(x, y, 'b', x, H(x), 'r')
#plt.show()
plt.figure()
plt.plot(t, a, 'bo', t, t-t, 'r')
plt.show()
plt.plot(t, a, 'bo', x, H(x), 'r')