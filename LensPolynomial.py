# -*- coding: utf-8 -*-
"""
Created on Thu May 21 12:40:31 2015

@author: Vladimir Putin
"""
import numpy as np
import matplotlib.pyplot as plt

"""
Here we are calculation polynomial coefficient describing
bent shape of a single capillary inside a Lens with
provided parameters
"""

# Make a picture explaining those
def getPolyCoeffs(x,y,D):
    y0 = y['y0']
    y1 = y['y1']
    ym = y['ym']
    y2 = y['y2']
    y3 = y['y3']
    h1 = x
    Din     = D['Din']
    Dout    = D['Dout']
    Dmax    = D['Dmax']
    h2 = h1 * Dout/Din
    hm = 0.5*Dmax * h1/(Din/2.)
    # y(x) = p0 + p1 x + p2 x**2 + p3 x**3 + p4 x**4 + p5 x**5
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

# NOTE - name __main__ is used only if the file is executed as itself
# not when it's imported
if __name__ == '__main__':
    y = {'y0' : 0., 'y1' : 40., 'ym' : 85, 'y2' : 140, 'y3' : 155}
    D = {'Din' : 4.5, 'Dout' : 2.4, 'Dmax' : 8.0}
    x_in = -1.1
    p = getPolyCoeffs(x_in,y,D)
    print p
    x = np.linspace(y['y1'],y['y2'],1000)
    yfin = p[0] + p[1]*x + p[2]*x**2 + p[3]*x**3 + p[4]*x**4 + p[5]*x**5
    fig1 = plt.figure(1,figsize=(10,4))
    ax1 = plt.subplot(111, label='dupa')
    ax1.set_xlim([y['y0'],y['y3']])
    ax1.set_ylim([-4, 4])
    ax1.plot(x,yfin,'r-', lw=2)
    ax1.plot([y['y0'],y['y1']],[0,x_in],'k-',lw=0.5)
    x2 = x_in * D['Dout'] / D['Din']
    ax1.plot([y['y2'], y['y3']],[x2,0],'k-',lw=0.5)
    plt.show()
