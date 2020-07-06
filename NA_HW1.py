# -*- coding: utf-8 -*-
"""
Created on Fri Aug 31 11:49:13 2018

@author: mechd
"""

def NA_HW1(n):
    import numpy as np
    import matplotlib.pyplot as plt
    k = np.empty(n+1)
    k_b = np.empty(n+1)
    N = np.arange(n+1)
    k_b[n] = 0
    k[0] = np.log(6) - np.log(5)
    for i in range(n):
        k[i+1] = -5*k[i] + 1/(i+1)
        k_b[n-i-1] = 0.2*(1/(n-i) - k_b[n-i])
        print('x' + str(i) + '=', k[i])
#        print('x_b' + str(n-i) + '=', k_b[n-i])
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.plot(k)
    ax.set_title('Forward iteration', fontsize=20)
    ax.set_xlabel('Iterations', fontsize=18)
    ax.set_ylabel('x_n', fontsize=16)
    for i,j in zip(N,k):
        ax.annotate(str(j),xy=(i,j))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.plot(k_b)
    ax.set_title('Backward iteration', fontsize=20)
    ax.set_xlabel('Iterations', fontsize=18)
    ax.set_ylabel('x_n', fontsize=16)
    for i,j in zip(N,k_b):
        ax.annotate(str(j),xy=(i,j))