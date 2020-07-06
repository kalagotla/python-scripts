# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 10:27:05 2017

@author: mechd
"""

def Amit2(i_n, jn, tn):
    import numpy as np
    delx = 0.01
    dely = 0.002
    delt = 0.1
    Tx1 = 200
    Tx0 = 400
    d_phi_d_x = (Tx1-Tx0)/(delx*i_n)
    Ty1 = 500
    Ty0 = 0
    d_phi_d_y = (Ty1-Ty0)/(dely*jn)
    phi_ij1 = 10
    hi_ij11 = 0
    phi_i11j = 0
    for t in range(tn):
        for j in range(jn):
            for i in range(i_n):
                phi_i01j = phi_ij1 - delx*d_phi_d_x
                phi_ij01 = phi_ij1 - dely*d_phi_d_y
                phi_ij2 = ((phi_i01j - 2*phi_ij1 + phi_i11j)/delx**2
                            + (phi_ij01 - 2*phi_ij1 + phi_ij11)/dely**2)*delt + phi_ij1
                phi_ij1 = phi_ij2
                
                print(phi_ij2)
    return