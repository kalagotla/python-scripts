# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 16:33:01 2018

@author: mechd
"""

def data_filter(infile, xmin, xmax, zmin, zmax):
    import numpy as np
    import math
    import matplotlib.pyplot as plt
    import time
    import scipy.interpolate
#    import plotly as py
#    import plotly.graph_objs as go
    start = time.time()
    
    #Read file into an array
    f_in = np.loadtxt(infile, skiprows=1)
    rows = len(f_in)
    columns = len(f_in[0])
    print('No. of rows:', rows)
    print('No. of columns:', columns)
    f_out = np.zeros((rows, columns))
    k = 0
    for i in range(rows):
        if xmin<=f_in[i, 0]<=xmax and zmin<=f_in[i,1]<=zmax:
            for j in range(columns):
                f_out[k,j] = f_in[i,j]
            k += 1
    k = k-1
    print('No. of rows in filtered plane:', k)
    
    np.savetxt('Datafile.txt', f_out[~np.all(f_out == 0, axis=1)])
    stop = time.time()
    
    xi, zi = np.linspace(xmin, xmax, 1000), np.linspace(zmin, zmax, 1000)
    xi, zi = np.meshgrid(xi, zi)
        
    # Interpolate; there's also method='cubic' for 2-D data such as here
    mi = scipy.interpolate.griddata((f_out[:,0], f_out[:,1]), f_out[:,2], (xi, zi), method='cubic')
    plt.figure()
#    plt.contourf(xi, zi, mi, cmap='jet')
    plt.pcolor(xi, zi, mi, cmap='jet')
    print('Elapsed time:', stop-start)
    
    return