# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 14:16:21 2016

@author: mechd
"""

def isNarcissistic(x):
    p = 0
    NumofDigits = len(str(x))
    for i in str(x):
        n = int(i)**NumofDigits
        p = p + n
    if p == x:
        return(True)
    else:
        return(False)
        
        
    