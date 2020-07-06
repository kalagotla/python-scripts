# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 01:20:46 2016

@author: mechd
"""

def SWPI1(M1, rhof1, T1, rhop, dp, R, gamma, mu, cd, n):
    import numpy as np
    import scipy.integrate as intg
    import matplotlib.pyplot as plt
    from scipy.interpolate import interp1d
    from scipy.interpolate import InterpolatedUnivariateSpline
    Vp = (4*np.pi*(dp/2)**3)/3
    P1 = rhof1*R*T1
    P2 = P1*(2*gamma*M1**2/(gamma+1) - ((gamma-1)/(gamma+1)))
    T2 = T1*((1 + 0.5*(gamma-1)*M1**2)*(2*gamma*M1**2/(gamma-1) - 1)*2*(gamma-1)/((gamma+1)**2*M1**2))
    rhof2 = rhof1*(P2*T1)/(P1*T2)
    vf1 = M1*np.sqrt(gamma*R*T1)
    vf2 = vf1*rhof1/rhof2
    M2 = vf2/np.sqrt(gamma*R*T2)
    # return vf1, vf2
    "SW eqns"
    Pr = 0.75
    G = (gamma+1)/((4/3)+((gamma-1)/Pr))
    X_SW = (8/G)*(mu/(vf1-vf2)*rhof1)
    X_SW_array = np.linspace(-X_SW, X_SW, n)
    # return X_SW_array
    f = (8/3)*(rhof1*vf1/mu)*((gamma-1)/gamma)*X_SW_array
    vf_SW = (vf1 + vf2*np.e**f) / (1 + np.e**f)
    t_SW = X_SW/vf_SW[5]
    # return vf_SW
    t = np.linspace(0, t_SW, n)
    # vf_SW_inter = interp1d(t, vf_SW, kind='cubic', fill_value='extrapolate')
    vf_SW_inter = InterpolatedUnivariateSpline(t, vf_SW, k=3)
    vf_SW_diff = vf_SW_inter.derivative()
    # return vf_SW_diff(5)
    # return vf_SW_diff(t)
    # return vf_SW
    def Fd(vp,t):
        # dvpdt = -3*np.pi*mu*dp*(vp[0]-vf_SW_inter)/(Vp*rhop)
        # dvpdt = -3*np.pi*mu*dp*(vp[0]-vf_SW_inter(t))/(Vp*rhop)
        dvpdt = (0.5*rhof1*0.25*np.pi*(dp**2)*cd*(vp[0]-vf_SW_inter(t))**2)/(Vp*rhop)
        return dvpdt
    def Fp(vp, t):
        dvpdt = (rhof1/rhop)*(vf_SW_diff(t))
        return dvpdt
    def Fam(vp, t):
        dvpdt = (0.5*rhof1*vf_SW_diff(t))/(rhop+0.5*rhof1)
        return dvpdt

    vp0 = vf1
    vp_d = intg.odeint(Fd, vp0, t)
    vp_p = intg.odeint(Fp, vp0, t)
    vp_am = intg.odeint(Fam, vp0, t)
    plt.plot(X_SW_array, Fd(vp_d, t), '-r', lw=5.0)
    plt.plot(X_SW_array, Fp(vp_p, t), '-b', lw=5.0)
    plt.plot(X_SW_array, Fam(vp_am, t), '-g', lw=5.0)
    plt.show()
    plt.figure()
    # return Fd(vp_d, t), vp_d, vf_SW
    # return Fd(vp_p, t), vp_p, vf_SW
    #return Fam(vp_am, t), Fp(vp_p, t), Fd(vp_d, t)
