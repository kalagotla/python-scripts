# -*- coding: utf-8 -*-
"""
Created on Wed Mar 28 08:22:40 2018

@author: mechd
"""
def PIV(filename, tol):
    import numpy as np
    import scipy.ndimage
    import matplotlib.pyplot as plt
    
    def isfloat(value):
        try:
            float(value)
            return True
        except ValueError:
            return False
    
    def Temp_files(filename):
        f = open(filename, 'r')
        g = open('Temp1.dat', 'w')
        h = open('Temp2.dat', 'w')
        i, j = 0, 0
        for line in f:
            lines = line
            columns = line.split()
            if i < 172*130 and isfloat(columns[1]) == True:
                g.write(str(i))
                g.write(' ')
                g.write(str(lines))
                i += 1
            elif isfloat(columns[1]) == True:
                h.write(str(j))
                h.write(' ')
                h.write(str(lines))
                j += 1
        g.close()
        h.close()
        f.close()
        f1 = np.loadtxt('Temp1.dat')
        f2 = np.loadtxt('Temp2.dat')
        return f1, f2
    
#    def Contours(filename, u_norm):
#        from scipy.ndimage.filters import gaussian_filter
#        f1, f2 = Temp_files(filename)
#        f1 = gaussian_filter(f1, 2.0)
##        f1 = scipy.ndimage.zoom(f1, 3)
#        plt.figure(1)
#        plt.contourf(f1, cmap='jet')
##        plt.imshow(f1, interpolation='gaussian', cmap='jet')
#        plt.colorbar()
#        return
#    
#    Contours(filename, u_norm)


    def Combine_data(filename, tol):
        f1, f2 = Temp_files(filename)
#        print(f1[171,1], f1[2*171,1])
        imax, jmax = 171, 129
        i = 0
        p1, p2 = 5,5
        error = tol*f1[:,3].max()
        v1_save, v2_save = 1000, 0
        print('error=', error)
#        delta_x, delta_y = np.zeros(100), np.zeros(100)
        v1, v2 = np.zeros(1000**2), np.zeros(1000**2)
        for i1 in range(50,80):
            for i2 in range(45,85):
                for i3 in range(20):
                    if (abs(f1[(i1+1)*(imax+1)-1, 3] - f2[i2*(imax+1)+i3, 3]) <= error 
                        and abs(f1[(i1+1)*(imax+1)-1, 3] - f2[i2*(imax+1)+i3, 3])!=0
#                        Check p1 addition to i1 & i2. lly with p2
                        and abs(f1[(i1+p1)*(imax+1)-1, 3] - f2[(i2+p1)*(imax+1)+i3, 3]) <= error
                        and abs(f1[(i1-p2)*(imax+1)-1, 3] - f2[(i2-p2)*(imax+1)+i3, 3]) <= error):
#                        print(i, i1, i2, i3,
#                              f1[(i1+1)*(imax+1)-1,3], f2[i2*(imax+1)+i3, 3])
                        v1[i] = f1[(i1+1)*(imax+1)-1,3]
                        v2[i] = f2[i2*(imax+1)+i3, 3]
                        i += 1
                        if (abs(v1_save-v2_save) >= abs((v1[i-1]-v2[i-1]))):
                            v1_save = v1[i-1]
                            v2_save = v2[i-1]
                            delta_x = f2[i2*(imax+1)+i3,1] - f1[(i1+1)*(imax+1)-1,1]
                            delta_y = f2[i2*(imax+1)+i3,2] - f1[(i1+1)*(imax+1)-1,2]
                            print('------------Saved above value---------------')
#                        break
#        k = 8
        f2[:,1] = f2[:,1] - delta_x
        f2[:,2] = f2[:,2] - delta_y
        f = open('AvgV1.dat', 'w')
        f.write('TITLE = "AvgVx"\nVARIABLES = "x", "y", "Avg Vx"\nZONE T="Frame 1" I=172, J=130, F=POINT\n')
        f.close()
        f = open('AvgV1.dat', 'ab')
        np.savetxt(f, f1[:,1:4], fmt='%f')
        f.close()
        f = open('AvgV1.dat', 'a')
        f.write('ZONE T="Frame 2" I=172, J=130, F=POINT\n')
#        f = open('AvgV2.dat', 'w')
#        f.write('TITLE = "AvgVx"\nVARIABLES = "x", "y", "Avg Vx"\nZONE T="Frame 1" I=172, J=130, F=POINT\n')
        f.close()
        f = open('AvgV1.dat', 'ab')
        np.savetxt(f, f2[:,1:4], fmt='%f')
        f.close()
        return
    
    Combine_data(filename, tol)
    
    return

PIV('AvgVx_2.DAT', 0.1)