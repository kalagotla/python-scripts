# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 10:29:27 2017

@author: mechd
"""

#n, dp = statistical_data(1000, 1, 1)

def statistical_data(total_par, del_dp, guess):
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.optimize import fsolve
    mu = 265
    sigma_g = 1.68
    mu_g = 229
    dp_min, dp_max = 213, 293
    dp_n = (dp_max-dp_min)/del_dp + 1
    dp_n = int(dp_n)
    di = np.linspace(213, 293, dp_n)
    def function(sigma):
        func1 = 0
        for i in range(dp_n):
            func1 += np.exp(-(di[i]-mu)**2 / (2*sigma**2)) / np.sqrt(2*np.pi*sigma**2)
            #func1 += (np.exp(-((di[i]-mu)**2)/(2*sigma**2))) * np.log(di[i])
            #func1 += (np.exp(-((di[i]-mu)**2)/(2*sigma**2)))*(np.log((di[i])/mu_g))**2
        #func2 = -np.sqrt(2*np.pi*sigma**2)*np.log(mu_g)
        #func2 = - np.sqrt(2*np.pi*sigma**2)*(np.log(sigma_g))**2
        #func = func1+func2
        func = func1-1
        return func
    sigma = np.linspace(0.1, 500, 10000)
    plt.figure(1)
    plt.plot(sigma, function(sigma))
    plt.show()
    print('\nsigma=', sigma, '\nfunc=', function(sigma))
    
    sigma_initial_guess = guess
    sigma_solution = fsolve(function, sigma_initial_guess)
    n = np.zeros(dp_n)
    dp = np.zeros(dp_n+1)
    dp[0] = dp_min
    for i in range(dp_n):
        n[i] = total_par * (np.exp(-(di[i]-mu)**2 / (2*sigma_solution**2))
                / np.sqrt(2*np.pi*sigma_solution**2))
        dp[i+1] = dp[i] + del_dp
    print('n=', n)
    print('sigma_solution=', sigma_solution)
    return n, dp
        