# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 15:47:22 2016

@author: mechd
"""

def FindNarcissistic(n):
    x = []
    for m in range(1, 10**9):
        if len(x) <= n:
            p = 0
            NumofDigits = len(str(m))
            for i in str(m):
                i = int(i)**NumofDigits
                p = p + i
            if p == m:
                x += list([m])
            if len(x) == n:
                return x