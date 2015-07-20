# -*- coding: utf-8 -*-
"""
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from PlotMono import plot2D
import screening as scr
from simulate_spoly import HexStructure
import multiprocessing as mp

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

# ray traycing settings (powerful pc defaults)    
mGlass  = rm.Material(('Si', 'O'), quantities=(1, 2), rho=2.2)
mGold   = rm.Material('Au', rho=19.3)
repeats = 5e4           # number of ray traycing iterations
processes = 8           # number of processes used
E0      = 9000.         # energy in electronoVolts
nRefl   = 125           # number of reflections

# Lower expectations for home computers
if mp.cpu_count() <= 2:
    repeats = 100
    processes = 1
    print "Running on a slow machine"
else:
    print "Running on a fast machine"

# Constant capillary/setup parameters [mm]
y0 =    0.      # relative light source position
y1 =    40.     # lens entrance
y2 =    140.    # lens exit
yf =    155.    # focus position
ym =    88.     # capillaries turning point
hMax =  4.0     # maximum possible distance from y = 0 axis
Din =   4.5     # lens entrance diameter
Dout =  2.4     # lens exit diameter
Dmax =  2*hMax  # max diameter
rIn =   0.006*30     # capillary radius
rOut = Dout/Din * rIn # Radius must shrink alongside the lens
rMax = Dmax/Din * rIn # Max value of local radius
wall =   0.001 # |*50 make wider walls for structure visibility

# Hex structure parameters
nx_capillary = 3
ny_bundle = 3

# Pinhole parameters
pinlen  = 0.005                 # Length 
rpin    = rIn / 2.0             # Pinhole radius [mm]
ypin    = 155.0 - pinlen        # Optical path position

# Source parameters
distx       = 'flat'
dx          = 0.1
distxprime  = 'flat'
dxprime     = 0.1
# z-direction
distz       = 'flat'
dz          = 0.1
distzprime  = 'flat'
dzprime     = 0.1

class BentCapillary(roe.OE):
    def __init__(self, *args, **kwargs):
        self.p  = kwargs.pop("curveCoeffs")
        self.rIn    = kwargs.pop("rIn")
        self.rOut   = kwargs.pop("rOut")
        self.rMax   = kwargs.pop("rMax")
        roe.OE.__init__(self, *args, **kwargs)
        self.y1 = self.limPhysY[0]
        self.y2 = self.limPhysY[1]
        self.ym = self.bl.ym

        # local radius function parameters (quadratic from rIn ot rOut)
        # linear solve of 3 equations for parabolic shape
        # of r(s): currently this curve is the same for each 
        # capillary, so this could've been done once for lens,
        # not once for capillary
        B = [self.rIn, self.rOut, self.rMax]
        A = []
        A.append([1., y1, y1**2])
        A.append([1., y2, y2**2])
        A.append([1., ym, ym**2])
        # Polynomial parameters for local radius-position relation
        self.pr = np.linalg.solve(A,B)
        self.isParametric = True

    def local_x0(self, s):
        return self.p[0] + self.p[1]*s + self.p[2]*s**2 + self.p[3]*s**3 + self.p[4]*s**4 + self.p[5]*s**5

    def local_x0Prime(self, s):
        return self.p[1] + 2*self.p[2]*s + 3*self.p[3]*s**2 + 4*self.p[4]*s**3 + 5*self.p[5]*s**4

    def local_r0(self, s):
        return self.pr[0] + self.pr[1]*s + self.pr[2]*s**2

    def local_r0Prime(self,s):
        return self.pr[1] + 2*self.pr[2]*s

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
    beamLine.ym = ym
    beamLine.hMax   = hMax
    beamLine.h1     = Din/2
    beamLine.Dout   = Dout
    beamLine.nRefl  = nRefl

    # [0] - Source of light
    rs.GeometricSource(
        beamLine,'GeometricSource',(0,0,0), nrays=nrays,
        distx=distx, dx=dx, distxprime=distxprime, dxprime=dxprime,
        distz=distz, dz=dz, distzprime=distzprime, dzprime=dzprime,
        distE='lines', energies=(E0,), polarization='horizontal')

    # Insert screen at the lens entrance here
    beamLine.entScreen = rsc.Screen(beamLine, 'EntranceScreen',(0,y1,0))

    # [1] - Lens related containers
    beamLine.capillaries = []
    beamLine.toPlot = []

    # Create object with (x,y) points describing hexagonal 
    # structure of polycapillary optics
    entrance_Structure = HexStructure(nx_capillary=nx_capillary, \
                                    ny_bundle=ny_bundle, \
                                    capillary_diameter=2*(rIn + wall))

    # Show obtained structure and save as png
    entrance_Structure.test()

    # Iterate through polar coordinates of those capillaries 
    # provided by a pythonic (?) generator
    for r, phi in entrance_Structure.genPolars():
        roll = phi
        x = r

        p = getPolyCoeffs(y0,y1,ym,y2,yf,x,Din,Dout,hMax)
        capillary = BentCapillary(beamLine, 'BentCapillary',
                [0,0,0], roll=roll, limPhysY=[y1, y2], order=8,
                rIn=rIn, rOut=rOut, rMax=rMax, material=mGlass,
                curveCoeffs=p)
        beamLine.capillaries.append(capillary)

        # Save capillaries shown on z=0 coss-section ? 
        # Z = 0 is no longer special
        # and as it is clear neither is phi = pi/3,
        # so some other idea for crosssection plot
        # is needed TODO
        # DEBUG quick polar to cartesian re-transformation
        x_cap = r * np.cos(phi)
        if abs(x_cap) < 0.005:
            beamLine.toPlot.append(len(beamLine.capillaries))

    print 'Number of capillaries: ' + str(len(beamLine.capillaries))

    beamLine.exitScreen = rsc.Screen(beamLine,'ExitScreen', (0,y2,0))

    # Pinhole structure related container(s)
    beamLine.pinholes = []
    
    # Insert very short and very thin golden and lined capillaries 
    # into the focus, acting as a proper one way image sharpening pinhole

    # FIXME - this is wrong!
    # Straight-Bent capillary polynomial factors:
    p_pin = [0, 0, 0, 0, 0, 0]

    # Focus size radius estimation
    focus_r = 0.09  # ? 

    for it in np.linspace(-1,1,11):
        x_in = it * focus_r
        # FIXME - this has to be changed into proper straight and short capillary
        # Because as it turns out h_in parameter is not producing expected behavior
        pinhole = BentCapillary(beamLine, 'PinHole',
                        [0,0,0], roll=0, limPhysY=[ypin, ypin+pinlen],
                        order=8, rIn=rpin, rOut=rpin, rMax=rpin,
                        material=mGold,
                        curveCoeffs=p_pin)
        beamLine.pinholes.append(pinhole)

    # Helpful print
    print 'Number of pinholes: ' + str(len(beamLine.pinholes))

    # Create evenly distributed screens between lens exit
    # and M=1 spot
    scr.createScreens(beamLine,[y2, yf + yf-y2], 3)
    # Set screen used before the pinhole
    scr.setUsed(beamLine, [0, ypin])
    # Set used after pinhole as well
    scr.setUsed(beamLine, [ypin, 200])

    return beamLine

def run_process(beamLine, shineOnly1stSource=False):
    # [0]
    beamSource = beamLine.sources[0].shine()
    # at the entrance | unused
    EntranceScreen = beamLine.entScreen.expose(beamSource)
    outDict = {'beamSource': beamSource, 'EntranceScreen': EntranceScreen}

    # [1]
    # Start collecting capillaries' light
    beamCapillaryGlobalTotal = None
    for i, capillary in enumerate(beamLine.capillaries):
        # Get both types of coordinates (global,local)
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
    prePinhole = scr.exposeScreens(beamLine, beamCapillaryGlobalTotal,\
            [0, ypin])

    # [2] - Light hits pinholes
    pinholeGlobalTotal = None
    for i, pinhole in enumerate(beamLine.pinholes):
        pinholeGlobal, pinholeLocal = pinhole.multiple_reflect(\
                beamCapillaryGlobalTotal, maxReflections=5)

        if pinholeGlobalTotal is None:
            pinholeGlobalTotal = pinholeGlobal
        else:
            good = (( pinholeGlobal.state == 1 ) |
                    ( pinholeGlobal.state == 3))
            # Accumulate photons at global beam
            rs.copy_beam(pinholeGlobalTotal, pinholeGlobal, good, includeState=True)

    # Expose screens to post pinholes beam
    postPinhole = scr.exposeScreens(beamLine, pinholeGlobalTotal,\
            [ypin, 200])

    # [2] - bypass the pinhole
    postNoPinhole = scr.exposeScreens(beamLine, beamCapillaryGlobalTotal,\
            [ypin, 200])

    # 
    outDict.update(prePinhole)
    outDict.update(postPinhole)
#    outDict.update(postNoPinhole)

    return outDict

rr.run_process = run_process

def main():
    beamLine = build_beamline()
    plot2D(beamLine)

    # Create xrtp.Plots in outside module  
    plots = scr.createPlots(beamLine, save=True)
    # FIXME: Manually add plot showing entrance structure
    limits = [ - Din/1.9, Din/1.9 ]
    plot = xrtp.XYCPlot(
        'EntranceScreen', (1, 3,),
        xaxis=xrtp.XYCAxis(r'$x$', 'mm', bins=256, ppb=2, limits=limits),
        yaxis=xrtp.XYCAxis(r'$z$', 'mm', bins=256, ppb=2, limits=limits),
        caxis='category', beamState='Screen_at_'+str(int(y2)))
    plot.title = 'Entrance structure'
    plot.baseName = 'Entrance_structure'
    plot.saveName = 'png/' + plot.baseName + '.png'
    plot.persistentName = 'pickle/' + plot.baseName + '.pickle'
    plots.append(plot)
    xrtr.run_ray_tracing(plots, repeats=repeats, beamLine=beamLine,\
            processes=processes)

    print 'Number of used capillaries: ', len(beamLine.capillaries)
#    return beamLine

if __name__ == '__main__':
    main()
