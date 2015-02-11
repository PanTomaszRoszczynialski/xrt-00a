# -*- coding: utf-8 -*-
"""
Created on Fri Feb 06 16:26:54 2015

@author: Vladimir Putin
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import xrt.backends.raycing as raycing
import xrt.backends.raycing.sources as rs
#import xrt.backends.raycing.apertures as ra
import xrt.backends.raycing.oes as roe
import xrt.backends.raycing.run as rr
import xrt.backends.raycing.materials as rm
import xrt.plotter as xrtp
import xrt.runner as xrtr
import xrt.backends.raycing.screens as rsc

mGlass = rm.Material(('Si', 'O'), quantities=(1, 2), rho=2.2)

E0 = 9000.
rSample = 100. # starting position of the lens
f = 500. # y length in mm from foucs to the end of the lens
r0 = 0.0051
wall = 0.02
layers = 20 # number of hexagonal layers
nRefl = 12
nReflDisp = 12 # unused
xzPrimeMax = 3.

class BentCapillary(roe.OE):
    def __init__(self, *args, **kwargs):
        self.rSample = kwargs.pop('rSample')
        self.entranceAlpha = kwargs.pop('entranceAlpha')
        self.f = kwargs.pop('f') # times two maybe?
        self.r0in = kwargs.pop('rIn')
        self.r0out = kwargs.pop('rOut')
        roe.OE.__init__(self, *args, **kwargs)

        s0 = self.f - self.rSample * np.cos(self.entranceAlpha)
        self.a0 = -np.tan(self.entranceAlpha) / 2 / s0
        self.b0 = self.rSample * np.sin(self.entranceAlpha) - self.a0 * s0**2
        self.s0 = s0
        self.ar = (self.r0out-self.r0in) / s0**2
        self.br = self.r0in
        self.isParametric = True

    def local_x0(self, s):  # axis of capillary, x(s)
        return self.a0 * (s)**2 + self.b0

    def local_x0Prime(self, s):
        return 2 * self.a0 * s

    def local_r0(self, s):  # radius of capillary (s)
        return self.ar * (s-self.s0)**2 + self.br

    def local_r0Prime(self, s):
        return self.ar * 2 * (s-self.s0)

    def local_r(self, s, phi):
        den = np.cos(np.arctan(self.local_x0Prime(s)))**2
        return self.local_r0(s) / (np.cos(phi)**2/den + np.sin(phi)**2)

    def local_n(self, s, phi):
        a = -np.sin(phi)
        b = -np.sin(phi)*self.local_x0Prime(s) - self.local_r0Prime(s)
        c = -np.cos(phi)
        norm = np.sqrt(a**2 + b**2 + c**2)
        return a/norm, b/norm, c/norm

    def xyz_to_param(self, x, y, z):
        """ *s*, *r*, *phi* are cynindrc-like coordinates of the capillary.
        *s* is along y in inverse direction, started at the exit,
        *r* is measured from the capillary axis x0(s)
        *phi* is the polar angle measured from the z (vertical) direction."""
        s = self.f - y
        phi = np.arctan2(x - self.local_x0(s), z)
        r = np.sqrt((x-self.local_x0(s))**2 + z**2)
        return s, phi, r

    def param_to_xyz(self, s, phi, r):
        x = self.local_x0(s) + r*np.sin(phi)
        y = self.f - s
        z = r * np.cos(phi)
        return x, y, z


def build_beamline(nrays=1000):
    beamLine = raycing.BeamLine(height=0)
    rs.GeometricSource(
        beamLine, 'GeometricSource', (0, 0, 0), nrays=nrays,
        dx=0., dz=0., distxprime='annulus',
        distE='lines', energies=(E0,), polarization='horizontal')
    beamLine.fsm1 = rsc.Screen(beamLine, 'DiamondFSM1', (0, rSample, 0))
    beamLine.capillaries = []
    beamLine.firstInLayer = []
    beamLine.xzMax = 0
    for n in range(layers):
        if n > 0:
            ms = range(n)
            i6 = range(6)
        else:
            ms = 0,
            i6 = 0,
        beamLine.firstInLayer.append(len(beamLine.capillaries))
        for i in i6:
            for m in ms:
                x = 2 * (r0+wall) * (n**2 + m**2 - n*m)**0.5
                alpha = np.arcsin(x / rSample)
                roll1 = -np.arctan2(np.sqrt(3)*m, 2*n - m)
                roll = roll1 + i*np.pi/3.
                capillary = BentCapillary(
                    beamLine, 'BentCapillary', [0, 0, 0], roll=roll,
                    material=mGlass, limPhysY=[rSample*np.cos(alpha), f],
                    f=f, rSample=rSample, entranceAlpha=alpha, rIn=r0, rOut=r0)
                beamLine.capillaries.append(capillary)
                if beamLine.xzMax < capillary.b0:
                    beamLine.xzMax = capillary.b0
    print 'max divergence =', alpha
    beamLine.xzMax += 2 * r0
    print len(beamLine.capillaries)
    beamLine.sources[0].dxprime = 0, np.arcsin((2*n+1) * (r0+wall) / rSample)
#    beamLine.sources[0].dxprime = (np.arcsin((2*n-3) * (r0+wall) / rSample),
#        np.arcsin((2*n+1) * (r0+wall) / rSample))
#    beamLine.sources[0].dxprime = 0, np.arcsin(r0 / rSample)
    beamLine.fsm2 = rsc.Screen(beamLine, 'DiamondFSM2', (0, f, 0))
    return beamLine

def run_process(beamLine, shineOnly1stSource=False):
    beamSource = beamLine.sources[0].shine()
#    raycing.rotate_beam(
#        beamSource, yaw=-beamLine.capillaries[0].entranceAlpha)
    beamFSM1 = beamLine.fsm1.expose(beamSource)
    outDict = {'beamSource': beamSource, 'beamFSM1': beamFSM1}
    beamCapillaryGlobalTotal = None
    for i, capillary in enumerate(beamLine.capillaries):
        beamCapillaryGlobal, beamCapillaryLocalN =\
            capillary.multiple_reflect(beamSource, maxReflections=nRefl)
        beamCapillaryLocalN.phi /= np.pi
        if beamCapillaryGlobalTotal is None:
            beamCapillaryGlobalTotal = beamCapillaryGlobal
        else:
            good = ((beamCapillaryGlobal.state == 1) |
                    (beamCapillaryGlobal.state == 3))
            rs.copy_beam(beamCapillaryGlobalTotal, beamCapillaryGlobal,
                         good, includeState=True)
        outDict['beamCapillaryLocalN{0:02d}'.format(i)] = beamCapillaryLocalN
    outDict['beamCapillaryGlobalTotal'] = beamCapillaryGlobalTotal
    beamFSM2 = beamLine.fsm2.expose(beamCapillaryGlobalTotal)
    outDict['beamFSM2'] = beamFSM2
    return outDict
rr.run_process = run_process


def plot2D():
    beamLine = build_beamline()
    fig1 = plt.figure(1, figsize=(8, 6))
#    ax1 = plt.subplot(111, aspect='equal', label='1')
    ax1 = plt.subplot(111, aspect=50, label='1')
    ax1.set_title('Cross-section of polycapillary at $z$=0')
    ax1.set_xlabel(r'$y$ (mm)', fontsize=14)
    ax1.set_ylabel(r'$x$ (mm)', fontsize=14)
    seq = [2, 5, 12, 5]
    for i in beamLine.firstInLayer:
        print i
        capillary = beamLine.capillaries[i]
        s = np.linspace(0, capillary.s0, 200)
        x = capillary.local_x0(s)
        r = capillary.local_r0(s)
        ax1.plot([0, f-s[-1]], [0, x[-1]], 'k-', lw=0.5)
        line = ax1.plot(f-s, x, 'k-.', lw=0.5)
        line[0].set_dashes(seq)
        ax1.plot(f-s, x-r, 'r-', lw=2)
        ax1.plot(f-s, x+r, 'r-', lw=2)
    ax1.set_xlim(0, f)
    ax1.set_ylim(-2*capillary.r0in, capillary.local_x0(0) + 2*capillary.r0out)
    fig1.savefig('PolycapillaryZ0crosssection.png')
    plt.show()
    
def main():
    beamLine = build_beamline()
    fwhmFormatStr3 = '%.3f'
    plots = []
#    PlotClass = xrtp.XYCPlotWithNumerOfReflections
    PlotClass = xrtp.XYCPlot

    limits = [-3,3]
    plot = xrtp.XYCPlot(
        'beamFSM2', (1, 3),
        xaxis=xrtp.XYCAxis(r'$x$', 'mm', bins=256, ppb=2, limits=limits),
        yaxis=xrtp.XYCAxis(r'$z$', 'mm', bins=256, ppb=2, limits=limits),
        caxis='category', beamState='beamFSM2', title='FSM1_Cat')
    plot.baseName = 'NCapillaries-a-FSM1Cat'
    plot.saveName = plot.baseName + '.png'
    plots.append(plot)
    xrtr.run_ray_tracing(plots, repeats=322, beamLine=beamLine,
                         processes=2)

#    for i in beamLine.firstInLayer:
#   capillariesToShow = 0, (layers-1)/2, layers-1    
    
if __name__ == '__main__':
    main()
#    plot2D()
