# -*- coding: utf-8 -*-
"""
Created on Wed Jul 25 14:08:34 2018

@author: mechd
"""

def vortex(n, K, r_min, r_max, z_range):
#    vortex(10, 733.038285, 0.01, 0.05, 100)
    # n defines fineness of the grid
    # K is the angular rotational velocity of the vortex.
    # It can be time dependent or constant. Here I'm using constant
    # r_min is minimum radius for polar coordinates. Below r_min there is no data
    # r_max is maximum radius for polar coordinates. Above r_max there is no data
    # z_range varies z from 0 to z_range
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
#    import plotly.plotly as py
#    import plotly.graph_objs as go
    x = np.linspace(-1, 1, n)
    y = np.linspace(-1, 1, n)
    X, Y = np.meshgrid(x, y)
    # Define all the variables after defining the grid
    r = (X**2 + Y**2)**0.5 # Center at the origin
    v = K/r**2 # Velocity definition for irrotational vortex
    rho = 1.273 # Density of the fluid
    p = -0.5*rho*v**2 # Pressure from bernoulli equation
    
    #Plot the data. This kind helps in plotting more plots in single figure.
    fig = plt.figure()
    ax = fig.add_subplot(211)
    contour = ax.contourf(X, Y, v, cmap='jet')
    ax.xaxis.grid(True)
    ax.yaxis.grid(True)
    ax.set_title('Velocity contours')
    ax1 = fig.add_subplot(212)
    contour1 = ax1.contourf(X, Y, p, cmap='jet')
    ax1.xaxis.grid(True)
    ax1.yaxis.grid(True)
    
    #Polar Co-ordinates
    r_range = np.linspace(r_min, r_max, n)
    theta = np.linspace(0, 2*np.pi, n)
    R, T = np.meshgrid(r_range, theta)
    v = K/R
    rho = 1.273 # Density of the fluid
    p = -0.5*rho*v**2 # Pressure from bernoulli equation
    # 2D polar contour
    fig = plt.figure()
    ax = fig.add_subplot(211, polar=True)
    contour = ax.contourf(T, R, v, cmap='jet')
    ax.set_title('Polar velocity contours')
    ax1 = fig.add_subplot(212, polar=True)
    contour1 = ax1.contourf(T, R, p, cmap='jet')
    
    # 3D contours
    z = np.linspace(-z_range, z_range, n)
    Z, Z1 = np.meshgrid(z, z)
    
    # fourth dimension - colormap
    # create colormap according to velocity
    color_dimension = v
    minn, maxx = color_dimension.min(), color_dimension.max()
    norm = matplotlib.colors.Normalize(minn, maxx)
    m = plt.cm.ScalarMappable(norm=norm, cmap='jet')
    m.set_array([])
    fcolors = m.to_rgba(color_dimension)
    
    # Plot on cartesian co-ordiantes
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    X, Y = R*np.cos(T), R*np.sin(T) # Grid transform onto cartesian co-ordiantes
    k1 = np.size(x)
    k2 = np.size(y)
    print(np.size(X), np.size(Y), k1, k2)
#    print(X, Y)
    f = np.empty((z_range*k1*k2,3)) #array to write into griddata
    g = np.empty((z_range*k1*k2,5)) #array to write into flowdata
    rho_data = np.full((k1, k2), rho) #Density data for flowdata
    v1 = -K*(Y/np.sqrt(X**2+Y**2)) #-u*sin(T)
    v2 = K*(X/np.sqrt(X**2+Y**2)) #u*cos(T)
#    v1 = K*(X**2 - Y**2) #X-component for velocity; Frank White example; pg no.: 263
#    v2 = -2*K*X*Y #Y-component for velocity; Same example
    v3 = np.zeros((k1, k2))
    v = np.sqrt(v1**2 + v2**2 + v3**2)
    p = -0.5*rho*v**2
    Temp = abs(p/(287.058*rho)) # Calculated from P = rho*R*T
    kc = 0
    kd = 0
    for i in range(z_range):
        Z = np.full((k1, k2), i)
#        Z, Z1 = np.meshgrid(z, z)
#        scatter = ax.scatter(X, Y, Z, cmap='jet')
#        surface = ax.plot_trisurf(x, y, z, cmap='jet')
        surface = ax.plot_surface(X, Y, Z, facecolors=fcolors)
        ax.set_title('Just grid with colors')
#        contour = ax.contourf3D(X, Y, Z, cmap='jet')
        for i in range(k1):
            for j in range(k2):
                f[kc, 0] = X[i, j]
                f[kc, 1] = Y[i, j]
                f[kc, 2] = Z[i, j]
                kc = kc + 1
                
        for i in range(k1):
            for j in range(k2):
                g[kd, 0] = rho_data[i,j]
                g[kd, 1] = v1[i, j]
                g[kd, 2] = v2[i, j]
                g[kd, 3] = v3[i, j]
                g[kd, 4] = Temp[i, j]
                kd = kd + 1
                
    fig = plt.figure()
    ax1 = fig.add_subplot(111, projection='3d')
    print(np.size(f), kc)
    ax1.scatter(f[0:z_range*k1*k2, 0], f[0:z_range*k1*k2, 1], f[0:z_range*k1*k2, 2], c='r')
    print('f_0_size=', np.size(f[0:z_range*k1*k2, 0]))
    ax1.set_xlabel('X Label')
    ax1.set_ylabel('Y Label')
    ax1.set_zlabel('Z Label')
        
    
    # Write data in a format to read into VISUAL3
    print(kc)
    np.savetxt('vortex_grid' + '.txt', f)
    np.savetxt('vortex_flow' + '.txt', g)
#    np.savetxt('test', Z)

    
    return