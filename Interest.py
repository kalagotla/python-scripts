# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 10:18:23 2017

@author: mechd
"""

def interest(mdy1, mdy2, principal, interest):
    import time
    import pdb
    import numpy as np
    import math
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import scipy.interpolate
    from mpl_toolkits.mplot3d import Axes3D
    import random
    from scipy.spatial import Delaunay
    from datetime import datetime
    #calculate number of days
    date_format = "%m/%d/%Y"
    a = datetime.strptime(mdy1, date_format)
    b = datetime.strptime(mdy2, date_format)
    delta = b - a
    n = delta.days+1
    print('Number of days = ', n)
    
    #calculate interest
    m = (n/30)*interest*principal
    print('Interest = ', m)
    
    
    