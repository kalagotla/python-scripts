# -*- coding: utf-8 -*-
"""
Created on Fri Jul 28 04:00:20 2007

@author: mechd
"""

def error_metric(xe, ye, Ue, Ve, xp, yp, Up, Vp, size):
    import numpy as np
    x_error = sum(sum(abs(xp[0:size, 0:size] - xe[0:size, 0:size])))/(size*size)
    y_error = sum(sum(abs(yp[0:size, 0:size] - ye[0:size, 0:size])))/(size*size)
    U_error = sum(sum(abs(Up[2:size, 0:size] - Ue[2:size, 0:size])))/((size-2)*size)
    V_error = sum(sum(abs(Vp[2:size, 0:size] - Ve[2:size, 0:size])))/((size-2)*size)
    
    print('x_error=', x_error, '\ny_error=', y_error, '\nU_error=', U_error, '\nV_error=', V_error)
    
    return