# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 18:52:51 2018

@author: mechd

for i in np.arange(0, 1, 0.05):
    f1, f2, g1, g2 = wave_equation(i, 100)
"""

def wave_equation(t, itr):
    import numpy as np
    import matplotlib.pyplot as plt
    import copy
    
    f1, g1 = np.zeros(2*itr), np.zeros(2*itr)
    f2, g2 = np.zeros(2*itr), np.zeros(2*itr)
#    u, f = np.zeros(10000), np.zeros(10000)
    def f_minus(t):
#        print('t=', t)
        j = 0
#        i = -1+t
#        if i<=-1+t:
#            f1[j] = 0
#            f2[j] = i
#            j = j+1
        for i in range(int((-1+t)*itr), int((0+t)*itr)):
            f1[j] = 1 + (i/itr-t)
            f2[j] = i/itr
#            print(i)
            j = j+1
        for i in range(int((0+t)*itr), int((0.5+t)*itr)):
            f1[j] = 1 - 2*(i/itr-t)
            f2[j] = i/itr
            j = j+1
#        i = 0.5+t
#        if i>=0.5+t:
#            f1[j] = 0
#            f2[j] = i
#            print(i)
#            j = j+1
#        return np.trim_zeros(f1, 'b'), np.trim_zeros(f2, 'b')
        return f1[0:j], f2[0:j]
    def f_plus(t):
        j = 0
#        i = -1-t
#        if i<=-1-t:
#            print('yes', j)
#            g1[j] = 0
#            g2[j] = i
#            j = j+1
        for i in range(int((-1-t)*itr), int((0-t)*itr)):
            g1[j] = 1 + (i/itr+t)
            g2[j] = i/itr
            j = j+1
        for i in range(int((0-t)*itr), int((0.5-t)*itr)):
            g1[j] = 1 - 2*(i/itr+t)
            g2[j] = i/itr
#            print(i)
            j = j+1
#        i = 0.5-t
#        if i>=0.5-t:
#            g1[j] = 0
#            g2[j] = i
#            j = j+1
#        return np.trim_zeros(g1, 'b'), np.trim_zeros(g2, 'b')
        return g1[0:j], g2[0:j]
    
#    t1 = np.linspace(-1, 0.5, np.size(f_minus(t)))
#    print('minus=', np.size(f_minus(t)), '\nplus=', np.size(f_plus(t)))


    f1, f2 = f_minus(t)
    g1, g2 = f_plus(t)
#    print(f1)
#    print(g1)
#    f_old, g_old = np.zeros(np.size(f1)), np.zeros(np.size(g1))
    f_old, g_old = copy.copy(f1), copy.copy(g1)
#    u = 0.5*np.append(f1, g1)
#    f = 0.5*np.append(f2, g2)
    for i in range(np.size(f2)):
        for j in range(np.size(g2)):
            if f2[i] == g2[j]:
                f1[i] = (f_old[i] + g_old[j])
#                print(0.5*f1[i], i, j)
                g1[j] = (f_old[i] + g_old[j])
#                print(0.5*g1[j], i, j)
                
#    u = np.trim_zeros(u, 'b')

#    print(0.5*f1)
#    print(0.5*g1)



    plt.plot(f2, 0.5*f1, 'b', g2, 0.5*g1, 'r', linewidth=5.0)
    plt.pause(0.5)
    plt.clf()
    
    return f1, f2, g1, g2