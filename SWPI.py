# -*- coding: utf-8 -*-
"""
Created on Sun Oct  9 12:10:04 2016

@author: mechd
"""
def SWPI(Mf, rhof, Tf, gamma, R, Vp, rhop, t, t2, Pgrad, Cd, mu, n):
    #SWPI(2.05, 1.273, 278, 1.4, 286.9, 3.053e-18, 1.3, 25, 50, 0, 0.4, 1.51e-5, 10)
    import numpy as np
    import scipy.integrate as intg
    import matplotlib.pyplot as plt
    from scipy.interpolate import interp1d
    Pf = rhof*R*Tf
    Pf2 = Pf*(2*gamma*Mf**2/(gamma+1) - ((gamma-1)/(gamma+1)))
    Tf2 = Tf*((1 + 0.5*(gamma-1)*Mf**2)*(2*gamma*Mf**2/(gamma-1) - 1)*2*(gamma-1)/((gamma+1)**2*Mf**2))
    rhof2 = rhof*(Pf2*Tf)/(Pf*Tf2)
    vf1 = Mf*np.sqrt(gamma*R*Tf)
    vf2 = vf1*rhof/rhof2
    Mf2 = vf2/np.sqrt(gamma*R*Tf2)
    t1 = np.linspace(0, t, n)
    t2 = np.linspace(t, t2, n)
    v_plot = vf1*np.ones(n)
    v2_plot = vf2*np.ones(n)
    "SW eqns"
    Pr = 0.75
    G = (gamma+1)/((4/3)+((gamma-1)/Pr))
    X_SW = (8/G)*(mu/(vf1-vf2)*rhof)
    # return X_SW
    X_SW_intg = np.linspace(0, X_SW, n)
    f = (8/3)*(rhof*vf1/mu)*((gamma-1)/gamma)*X_SW_intg
    vf_SW = (vf1 + vf2*np.e**f) / (1 + np.e**f)
    # vf_SW_inter = interp1d(t1, vf_SW, kind='cubic')
    t_SW = X_SW/vf_SW[n-7]
    # return t_SW
    "Particle eqns"
    rp = (3*Vp/(4*np.pi))**(1/3)
    # return rp
    def EOMP1(vp, t):
        Fp = Vp*(Pgrad)
        # Fam = 0.5*rhof*Vp*(vp-vf1)
        Fam = 0.5*rhof*Vp
        v = (vp-vf1)**2
        Fd = 2*np.pi*Cd*rhof*v*(rp**2)
        # Fd = 6*np.pi*mu*rp*(vp-vf1)
        F = rhop*Vp
        # return -(Fp+Fam+Fd)/F
        return -(Fp+Fd)/(F+Fam)
    def EOMP2(vp, t):
        Fp = Vp*(Pgrad)
        # Fam = 0.5*rhof*Vp*(vp-vf2)
        Fam = 0.5*rhof*Vp
        v = (vp-vf2)**2
        Fd = 2*np.pi*Cd*rhof*v*(rp**2)
        # Fd = 6*np.pi*mu*rp*(vp-vf2)
        F = rhop*Vp
        # return -(Fp+Fam+Fd)/F
        return -(Fp+Fd)/(F+Fam)
    vp_f = intg.odeint(EOMP1, vf1, t1)
    vp2_f = intg.odeint(EOMP2, vp_f[n-1], t2)
    # return vp_f
    "Path of particle"
    vp_inter = np.squeeze(vp_f)
    vp2_inter = np.squeeze(vp2_f)
    # return vp_inter
    vp_f_inter = interp1d(t1, vp_inter, kind='cubic')
    v_inter = interp1d(t1, v_plot, kind='cubic')
    vp2_f_inter = interp1d(t2, vp2_inter, kind='cubic')
    v2_inter = interp1d(t2, v2_plot, kind='cubic')
    # return vp_f_inter(5)
    x = np.zeros(n)
    xf = np.zeros(n)
    x2 = np.zeros(n)
    xf2 = np.zeros(n)
    # return t1[9]
    for i in range(n):
        x[i] = intg.quad(vp_f_inter, 0, t1[i])[0]
        xf[i] = intg.quad(v_inter, 0, t1[i])[0]
    # return x, xf
    for i in range(n):
        x2[i] = intg.quad(vp2_f_inter, t1[n-1], t2[i])[0]
        xf2[i] = intg.quad(v2_inter, t1[n-1], t2[i])[0]
        # x2[i] = x[n-1] + intg.quad(vp2_f_inter, t1[n-1], t2[i])[0]
        # xf2[i] = xf[n-1] + intg.quad(v2_inter, t1[n-1], t2[i])[0]
    plt.plot(t1, vp_f, t1, v_plot, t2, vp2_f, t2, v2_plot)
    plt.figure()
    plt.plot(t1, x, t1, xf, t2, x2, t2, xf2)
    plt.figure()
    plt.plot(t2, x2-xf2)
    plt.figure()
    plt.plot(xf2, vp2_f)
    print('Rp=', rp, '\nPf1=', Pf, '\nvf1=', vf1, '\nX_SW=', X_SW, '\nvf_SW=', vf_SW, '\nPf2=', Pf2, '\nTf2=', Tf2, '\nrhof2=', rhof2, '\nvf2=', vf2, '\nMf2=', Mf2)
    # return x, xf, x2, xf2
    return np.around(x-xf, 2), np.around(x2-xf2, 2)
    """x = np.zeros(10)
    xf = np.zeros(10)
    x2 = np.zeros(10)
    xf2 = np.zeros(10)
    for i in range(0, 10):
        x[i] = vp_f[i]*t1[i]
        xf[i] = v_plot[i]*t1[i]
        x2[i] = vp2_f[i]*t2[i]
        xf2[i] = v2_plot[i]*t2[i]
    plt.plot(t1, x, t1, xf, t2, x2, t2, xf2)
    print('vf1=', vf1, '\nPf2=', Pf2, '\nT2=', T2, '\nrho2=', rho2, '\nv2=', v2, '\nMf2=', Mf2)
    return vp_f, vp2_f, x-xf"""
    # return x, xf, x2, xf2, t1