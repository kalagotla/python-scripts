# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 12:54:47 2016

@author: mechd
"""
def Factorial(x):
    if x < 0:
        print("Error: Factorial is defined for non-negative integers only.")
    if x > 0:
        m=1
        for i in range(1,x+1):
            m = i * m
        return m
    if x == 0:
        return 1