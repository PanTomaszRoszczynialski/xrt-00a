# -*- coding: utf-8 -*-
"""
Created on Thu Apr 02 16:05:30 2015

@author: Vladimir Putin
"""
import numpy as np
#import matplotlib as mpl
import matplotlib.pyplot as plt

def plot2D(beamLine):
    fig1 = plt.figure(1,figsize=(8,6))
    ax1 = plt.subplot(111, aspect=20, label='1')
    ax1.set_title('Cross-section of polycapillary at z=0')
    ax1.set_xlabel('y [mm]', fontsize=14)
    ax1.set_ylabel('x [mm]', fontsize=14)
    seq = [2, 5, 12, 5]
    s = np.linspace(beamLine.y1, beamLine.y2,200)
    for capillary in beamLine.capillaries:
        x = capillary.local_x0(s)
        r = capillary.local_r0(s)
        ax1.plot([0, beamLine.y1],[0 , capillary.h_in],'k-', lw=0.5)
        ax1.plot([beamLine.y2, beamLine.yf],[x[-1] , 0],'k-', lw=0.5)
        ax1.plot(s, x+r,'r-', lw=1)
        ax1.plot(s, x-r,'r-', lw=1)
    ax1.set_xlim(0,beamLine.yf)
    ax1.set_ylim(-beamLine.hMax,beamLine.hMax)
    plt.show()
    fig1.savefig('Z0crosssection.png')
