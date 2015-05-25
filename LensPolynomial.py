# -*- coding: utf-8 -*-
"""
Created on Thu May 21 12:40:31 2015

@author: Vladimir Putin
"""
import numpy as np
from sympy import *
import matplotlib.pyplot as plt

#a = np.array([[3,1.3], [1,2]])
#b = np.array([9,8])
#x = np.linalg.solve(a, b)

# y(x) - lens describing function
# S1(x) - ray path from light source to lens entrance
# S2(x) - ray path from lens exit to the point of focus

# Make a picture explaining those
def getPolyCoeffs(y0,y1,ym,y2,y3,h1,Din,Dout,hMax):
    h2 = h1 * Dout/Din
    hm = hMax * h1/(Din/2.)
    # y(x) = p0 + p1 x + p2 x**2 + p3 x**3 + p4 x**$ + p5 x**5
    a = []
    b = []
    # Preparing A*p = B matrices for linalg.solve 
    # y(y1) == h1
    a.append(shape_coeffs(y1))
    b.append(h1)
    # y(y2) == h2
    a.append(shape_coeffs(y2))
    b.append(h2)
    # y(ym) == hm
    a.append(shape_coeffs(ym))
    b.append(hm)
    # y'(y1) == S1'(y1)
    # first find left slope
    s1 = h1/(y1-y0)
    a.append(diff_coeffs(y1))
    b.append(s1)
    # y'(y2) == S2'(x)
    s2 = -h2/(y3-y2)
    a.append(diff_coeffs(y2))
    b.append(s2)
    # y'(ym) == 0
    a.append(diff_coeffs(ym))
    b.append(0.)
    p = np.linalg.solve(a,b)
    return p

def shape_coeffs(x):
    coeff_array = [1., x, x**2, x**3, x**4, x**5]
    return coeff_array

def diff_coeffs(x):
    coeff_array = [0., 1., 2*x, 3*x**2, 4*x**3, 5*x**4]
    return coeff_array

if __name__ == '__main__':
    Din, Dout = 4.5, 2.4
    y0,y1,ym,y2,y3,h1 = 0., 40., 85., 140., 155., -1.1
    hMax = 4.
    h2 = h1 * Dout/Din
    p = getPolyCoeffs(y0,y1,ym,y2,y3,h1,Din,Dout,hMax)
    print p
    x = np.linspace(y1,y2,1000)
    y = p[0] + p[1]*x + p[2]*x**2 + p[3]*x**3 + p[4]*x**4 + p[5]*x**5
    fig1 = plt.figure(1,figsize=(10,4))
    ax1 = plt.subplot(111, label='dupa')
    ax1.set_xlim([y0,y3])
    ax1.set_ylim([-4, 4])
    ax1.plot(x,y,'r-', lw=2)
    ax1.plot([y0,y1],[0,h1],'k-',lw=0.5)
    ax1.plot([y2,y3],[h2,0],'k-',lw=0.5)
    plt.show()
