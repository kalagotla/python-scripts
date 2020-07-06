# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 15:22:27 2017

@author: mechd
"""

def OSWPI(Mf1, rhof1, Tf1, gamma, R, delta, n, tol):
#    OSWPI(1.5, 1.273, 300, 1.4, 286.9, 0.161443, 1000000, 1e-3)
#    OSWPI(2.75, 0.123, 114.317, 1.4, 287, 0.135263, 10000000, 1e-3)
    import numpy as np
    import matplotlib.pyplot as plt
    for i in range(n):
        sigma = i*np.pi/n
        Pf1 = rhof1*R*Tf1
        Mfn1 = Mf1*np.sin(sigma)
        Pf2 = Pf1*(2*gamma*Mfn1**2/(gamma+1) - ((gamma-1)/(gamma+1)))
        Tf2 = Tf1*((1 + 0.5*(gamma-1)*Mfn1**2)*(2*gamma*Mfn1**2/(gamma-1) - 1)*2*(gamma-1)/((gamma+1)**2*Mfn1**2))
        rhof2 = rhof1*(Pf2*Tf1)/(Pf1*Tf2)
        vfn1 = Mfn1*np.sqrt(gamma*R*Tf1)
        vfn2 = vfn1*rhof1/rhof2
        vf1 = vfn1/np.sin(sigma)
        vf2 = vfn2/np.sin(sigma-delta)
        vft1 = vf1*np.cos(sigma)
        vft2 = vf2*np.cos(sigma-delta)
        Mfn2 = vfn2/np.sqrt(gamma*R*Tf2)
        Mf2 = Mfn2/np.sin(sigma-delta)
        Mft1 = Mf1*np.cos(sigma)
        Mft2 = Mf2*np.cos(sigma-delta)
        PT1 = Pf1/(gamma-1) + 0.5*vf1**2/rhof1
        PT2 = Pf2/(gamma-1) + 0.5*vf2**2/rhof2
        if abs(vft1-vft2) <= tol and Pf2 > Pf1:
#            return Pf1, Mfn1, Pf2, Tf2, rhof2, vfn1, vfn2, 
#            vf1, vf2, vft1, vft2, sigma, Mfn2, Mf2, Mft1, Mft2, PT1, PT2
            print('Pf1=', Pf1, '\nvf1=', vf1, '\nPf2=', Pf2, '\nTf2=', 
                  Tf2,'\nrhof1=', rhof1, '\nrhof2=', rhof2, '\nrhof1inv=', 
                  1/rhof1, '\nrhof2inv=', 1/rhof2, '\nvf2=', vf2, '\nMf2=', Mf2)
            print('vfn1=', vfn1, '\nvfn2=', vfn2, '\nMfn2=', Mfn2, '\nvft1=', 
                  vft1, '\nvft2=', vft2, '\nMft2=', Mft2, '\nMft1=', 
                  Mft1, '\nPT1=', PT1, '\nPT2=', PT2, '\nsigma=', sigma,
                  '\nsigma_degrees=', 180*sigma/np.pi, '\npi_sigma=', np.pi/sigma)
            print('-------------------------------------------------------')

    return

