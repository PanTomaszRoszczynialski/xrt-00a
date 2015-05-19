# -*- coding: utf-8 -*-
"""
Created on Thu Apr 02 16:05:30 2015

@author: Vladimir Putin
"""
import numpy as np
#import matplotlib as mpl
import matplotlib.pyplot as plt

def plot2D(beamLine,f):
#    beamLine = build_beamline()
    fig1 = plt.figure(1, figsize=(8, 6))
#    ax1 = plt.subplot(111, aspect='equal', label='1')
    ax1 = plt.subplot(111, aspect=50, label='1')
    ax1.set_title('Cross-section of polycapillary at $z$=0')
    ax1.set_xlabel(r'$y$ (mm)', fontsize=14)
    ax1.set_ylabel(r'$x$ (mm)', fontsize=14)
    seq = [2, 5, 12, 5]

    capillary = beamLine.capillary
    s = np.linspace(0, capillary.s0, 100)
    x = capillary.local_x0(s)
    r = capillary.local_r0(s)
    ax1.plot([0, f-s[-1]], [0, x[-1]], 'k-', lw=0.5)
    line = ax1.plot(f-s, x, 'k-.', lw=0.5)
    line[0].set_dashes(seq)
    ax1.plot(f-s, x-r, 'r-', lw=2)
    ax1.plot(f-s, x+r, 'r-', lw=2)
    ax1.set_xlim(0, f+20)
#    ax1.set_ylim(0.0,0.2)
    ax1.set_ylim(-2*capillary.r0in - 0.05, capillary.local_x0(50) + 2*capillary.r0out + 0.05)
    ax1.set_aspect('auto')
    ax1.relim()
    ax1.autoscale_view()
    fig1.savefig('PolycapillaryZ0crosssection.png')
    plt.show()
