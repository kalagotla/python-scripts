#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 10:58:25 2018

@author: kalagodk

statistical_data('Experimental_plane/Folw_data/flow_data_254.txt', 'Experimental_plane/Folw_data/flow_data_255.txt','Experimental_plane/Folw_data/flow_data_256.txt','Experimental_plane/Folw_data/flow_data_257.txt','Experimental_plane/Folw_data/flow_data_258.txt','Experimental_plane/Folw_data/flow_data_259.txt','Experimental_plane/Folw_data/flow_data_260.txt','Experimental_plane/Folw_data/flow_data_261.txt','Experimental_plane/Folw_data/flow_data_262.txt','Experimental_plane/Folw_data/flow_data_263.txt','Experimental_plane/Folw_data/flow_data_264.txt','Experimental_plane/Folw_data/flow_data_265.txt','Experimental_plane/Folw_data/flow_data_266.txt','Experimental_plane/Folw_data/flow_data_267.txt','Experimental_plane/Folw_data/flow_data_268.txt','Experimental_plane/Folw_data/flow_data_269.txt','Experimental_plane/Folw_data/flow_data_270.txt','Experimental_plane/Folw_data/flow_data_271.txt','Experimental_plane/Folw_data/flow_data_272.txt')
"""

def statistical_data(f1, f2, f3, f4, f5, f6, f7, f8, f9 ,f10, f11, f12, f13, f14, f15, f16, f17, f18, f19):
    import numpy as np
#    def simplecount(filename1, filename2):
#        n = 0
#        m = 0
#        for line in open(filename1):
#            n += 1
#        for line in open(filename2):
#            m += 1
#        n = n - 1
#        m = m - 1
#        return n, m
#    n, m = simplecount(filename1, filename2)
#    print('n=', n, 'm=', m)
    
#    def postscript(filename1, filename2, n, m):
    g1 = np.loadtxt(f1)
    g2 = np.loadtxt(f2)
    g3 = np.loadtxt(f3)
    g4 = np.loadtxt(f4)
    g5 = np.loadtxt(f5)
    g6 = np.loadtxt(f6)
    g7 = np.loadtxt(f7)
    g8 = np.loadtxt(f8)
    g9 = np.loadtxt(f9)
    g10 = np.loadtxt(f10)
    g11 = np.loadtxt(f11)
    g12 = np.loadtxt(f12)
    g13 = np.loadtxt(f13)
    g14 = np.loadtxt(f14)
    g15 = np.loadtxt(f15)
    g16 = np.loadtxt(f16)
    g17 = np.loadtxt(f17)
    g18 = np.loadtxt(f18)
    g19 = np.loadtxt(f19)
    
    h = ((2.59e-7)*(g1+g19) + (2.05e-5)*(g2+g18) + (9.68e-4)*(g3+g17) + (2.74e-2)*(g4+g16) +
         (4.63e-1)*(g5+g15) + (4.68)*(g6+g14) + (28.3)*(g7+g13) + (102.3)*(g8+g12)
         + (221.2)*(g9+g11) + (286.06)*(g10))/1000
    
    np.savetxt('average.txt', h)
    
    return