# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 22:36:04 2016

@author: mechd
"""
import math
import numpy as np 

for n in range(13, 16):
    for t in np.arange(1.5e-3, 3e-3, 0.5e-3):
        #t = 1.5e-3
        #n = 14
        Cd = 0.7
        r1 = 2.48e-2
        Asw = .01098
        theta = 60*2*math.pi/360
        spacing = ((2*math.pi*r1) - (n*t))/n
        H = Asw / (n*spacing*Cd)
        b = n*t/(2*r1*math.pi*math.cos(theta))
        S = r1*math.tan(theta)/(2*H*(1-b))
        print(S)