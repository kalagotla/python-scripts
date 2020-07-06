# -*- coding: utf-8 -*-
"""
Created on Fri Mar 23 16:31:18 2018

@author: mechd
"""

def Data_separator(filename, filename2, n_col):
    import numpy as np
    def isfloat(value):
        try:
            float(value)
            return True
        except ValueError:
            return False
#    n = 0
#    f = np.loadtxt(filename, skiprows=1)
    f = open(filename, 'r')
    g = open(filename2, 'w')
    for line in f:
#        n += 1
#        print(repr(line))
#        columns1 = line.split()
        columns = line
        
        for i in range(n_col):
            if isfloat(columns[i]) == True:
#                print(float(columns[i]))
                g.write(str(columns))
#            else:
#                print(columns.split())
#        g.write('\n')
#        print(columns)
#    print('n=', n)
#    print('f=', f)
    
    f.close()
    g.close()
    
    
    return