# -*- coding: utf-8 -*-
"""
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from PlotMono import plot2D
from screening import *

import xrt.backends.raycing as raycing
import xrt.backends.raycing.sources as rs
#import xrt.backends.raycing.apertures as ra
import xrt.backends.raycing.oes as roe
import xrt.backends.raycing.run as rr
import xrt.backends.raycing.materials as rm
import xrt.plotter as xrtp
import xrt.runner as xrtr
import xrt.backends.raycing.screens as rsc
from LensPolynomial import getPolyCoeffs

# ray traycing settings    
mGlass  = rm.Material(('Si', 'O'), quantities=(1, 2), rho=2.2)
repeats = 5e4           # number of ray traycing iterations
E0      = 9000.         # energy in electronoVolts
nRefl   = 125           # number of reflections

# Constant capillary/setup parameters [mm]
y0 =    0.      # relative light source position
y1 =    40.     # lens entrance
y2 =    140.    # lens exit
yf =    155.    # focus position
ym =    88.     # capillaries turning point
hMax =  4.0     # maximum possible distance from y = 0 axis
Din =   4.5     # lens entrance diameter
Dout =  2.4     # lens exit diameter
rIn =   0.006     # lens radius
wall=   0.0005

# Surce parameters
distx       = 'flat'
dx          = 0.1
distxprime  = 'normal'
dxprime     = 0.1
# z-direction
distz       = 'flat'
dz          = 0.1
distzprime  = 'normal'
dzprime     = 0.1

class BentCapillary(roe.OE):
    def __init__(self, *args, **kwargs):
        self.y2 = kwargs.pop("y2")
        self.p  = kwargs.pop("curveCoeffs")
        self.rIn    = kwargs.pop("rIn")
        self.h_in   = kwargs.pop("h_in")
        roe.OE.__init__(self, *args, **kwargs)
        self.isParametric = True

    def local_x0(self, s):
        return self.p[0] + self.p[1]*s + self.p[2]*s**2 + self.p[3]*s**3 + self.p[4]*s**4 + self.p[5]*s**5

    def local_x0Prime(self, s):
        return self.p[1] + 2*self.p[2]*s + 3*self.p[3]*s**2 + 4*self.p[4]*s**3 + 5*self.p[5]*s**4

    def local_r0(self, s):
        return self.rIn

    def local_r0Prime(self,s):
        return 0

    def local_r(self, s, phi):
        den = np.cos(np.arctan(self.local_x0Prime(s)))**2
        return self.local_r0(s) / (np.cos(phi)**2/den + np.sin(phi)**2)

    def local_n(self, s, phi):
        a = -np.sin(phi)
        b = -np.sin(phi)*self.local_x0Prime(s) - self.local_r0Prime(s)
        c = -np.cos(phi)
        norm = np.sqrt(a**2 + b**2 + c**2)
        # FIXME: a and c probably should also get minues
        # but due to symmetry it's not visible for now?
        return a/norm, -b/norm, c/norm

    def xyz_to_param(self, x, y, z):
        """ *s*, *r*, *phi* are cynindrc-like coordinates of the capillary.
        *s* is along y in inverse direction, started at the exit,
        *r* is measured from the capillary axis x0(s)
        *phi* is the polar angle measured from the z (vertical) direction."""
        s = y
        phi = np.arctan2(x - self.local_x0(s), z)
        r = np.sqrt((x-self.local_x0(s))**2 + z**2)
        return s, phi, r

    def param_to_xyz(self, s, phi, r):
        x = self.local_x0(s) + r*np.sin(phi)
        y = s
        z = r * np.cos(phi)
        return x, y, z

def build_beamline(nrays=1e4):
    beamLine = raycing.BeamLine(height=0)
    beamLine.y0 = y0
    beamLine.y1 = y1
    beamLine.y2 = y2
    beamLine.yf = yf
    beamLine.hMax   = hMax
    beamLine.h1     = Din/2
    beamLine.Dout   = Dout
    beamLine.nRefl  = nRefl

    rs.GeometricSource(
        beamLine,'GeometricSource',(0,0,0), nrays=nrays,
        distx=distx, dx=dx, distxprime=distxprime, dxprime=dxprime,
        distz=distz, dz=dz, distzprime=distzprime, dzprime=dzprime,
        distE='lines', energies=(E0,), polarization='horizontal')

    # Insert screen at the lens entrance here
    beamLine.entScreen = rsc.Screen(beamLine, 'EntranceScreen',(0,y1,0))

    beamLine.capillaries = []
    layers = range(10,20) + range(30,40) + range(90,100) + range(130,135)
    beamLine.toPlot = []
    for n in layers:
        if n > 0:
            ms = range(n)
            i6 = range(6)
        else: # this would happen only if layers were negative?
            ms = 0,
            i6 = 0,
        beamLine.toPlot.append(len(beamLine.capillaries))
        for i in i6:
            for m in ms:
                # this seems like h_in
                x = 2*(rIn + wall) * (n**2 + m**2 - n*m)**0.5
                roll1 = -np.arctan2(np.sqrt(3)*m, 2*n - m)
                roll = roll1 + i*np.pi/3.
                p = getPolyCoeffs(y0,y1,ym,y2,yf,x,Din,Dout,hMax)
                capillary = BentCapillary(beamLine, 'BentCapillary',
                        [0,0,0], roll=roll, limPhysY=[y1, y2], order=8,
                        rIn=rIn, curveCoeffs=p, y2=y2, h_in=x)
                beamLine.capillaries.append(capillary)
    print 'Number of capillaries: ' + str(len(beamLine.capillaries))

    beamLine.exitScreen = rsc.Screen(beamLine,'ExitScreen', (0,y2,0))

    # Create evenly distributed screens between lens exit
    # and M=1 spot
    createScreens(beamLine,[y2, yf + yf-y2], 11)
    return beamLine


def run_process(beamLine, shineOnly1stSource=False):
    beamSource = beamLine.sources[0].shine()
    # at the entrance | unused
    EntranceScreen = beamLine.entScreen.expose(beamSource)
    outDict = {'beamSource': beamSource, 'EntranceScreen': EntranceScreen}
    # Start collecting capillaries' light
    beamCapillaryGlobalTotal = None
    for i, capillary in enumerate(beamLine.capillaries):
        # Get both type of coordinates (global,local)
        beamCapillaryGlobal, beamCapillaryLocalN =\
            capillary.multiple_reflect(beamSource, maxReflections=nRefl)
        # Not sure what is this for 
        beamCapillaryLocalN.phi /= np.pi

        if beamCapillaryGlobalTotal is None:
            beamCapillaryGlobalTotal = beamCapillaryGlobal
        else:
            good = ((beamCapillaryGlobal.state == 1) |
                    (beamCapillaryGlobal.state == 3))
            # Add photons to GlobalTotal
            rs.copy_beam(beamCapillaryGlobalTotal, beamCapillaryGlobal,
                         good, includeState=True)

    # Prepare acces to Global beam 
    # (individual capillaries might be acessed as well)
    outDict['beamCapillaryGlobalTotal'] = beamCapillaryGlobalTotal

    # See them on screen 
    ExitScreen = beamLine.exitScreen.expose(beamCapillaryGlobalTotal)
    outDict['ExitScreen'] = ExitScreen
    # Create exposed beamlines in outside module screens.
    outDict.update(exposeScreens(beamLine, beamCapillaryGlobalTotal))

    return outDict

rr.run_process = run_process

def main():
    beamLine = build_beamline()
    plot2D(beamLine)

    # Create xrtp.Plots in outside module  
    plots = createPlots(beamLine)
    xrtr.run_ray_tracing(plots, repeats=repeats, beamLine=beamLine, processes=8, threads=8)

#    return beamLine

if __name__ == '__main__':
    main()
