# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 17:20:43 2016

@author: mechd
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 03:15:31 2016

@author: Kai
"""

import numpy as np
import math
from scipy.linalg import solve
import matplotlib.pyplot as plt
def fred(n):
  a=-1
  b=1
  lda=-1
  def g(x):
   return (np.exp(x))

  def K(x,y):
    return x*np.e**(y*(1-x))
    
  h=(b-a)/(n)
  
  c=b+h
  x1=np.arange(a,c,h)
  B=g(x1)
  K1 = np.zeros((n+1, n+1))
  for i in range(n+1):
            for j in range(n+1):
                K1[i, j] = K(a+h*i, a+h*j)
  w=[]
  w.append(h/2)
  i=2
  while i<=n:
      w.append(h)
      i=i+1
  w.append(h)
  
  W= np.diag(w)
  
  i=np.ones(n+1)
  I=np.diag(i)
  A=I-(np.dot(K1,W))
  Ain=np.linalg.inv(A)
  z=solve(A,B)
  plt.plot(x1,z,'k')
  plt.show()
  return z,x1

def error(n1,n2):
  u=0  
  q1,o1=fred(n1)
  q2,o2=fred(n2)
  print (len(q2))
  e=-1
  
  while u<=n1:
      #if abs(q1[u]-((q2[u]+q2[u+1])/2))>=e:   
      e=abs(q1[u]-((q2[u]+q2[u+1])/2))
      print (e)
      
      u=u+1  
  return e
def convergenceratio():
    del12=error(100,200)
    del24=error(200,400)
    del48=error(400,800)
    del816=error(800,1600)
    r1=del12/del24
    r2=del24/del48
    r3=del48/del816
    
    return r1,r2,r3
    