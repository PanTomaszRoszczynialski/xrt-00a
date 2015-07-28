# -*- coding: utf-8 -*-
"""
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import glob
import os

from CapillaryElements import PolyCapillaryLens, StraightCapillary
from CapillaryElements import Pinhole

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
repeats = 5e4           # number of ray traycing iterations
processes = 8           # number of processes used
E0      = 9000.         # energy in electronoVolts
nRefl   = 125           # number of reflections
save    = True          # save results as pickles?

# Delete all pickle files (they can always be recovered from git)
picklePaths = glob.glob('pickle/*.pickle')

for path in picklePaths:
    print 'Deleting pickle:', path
    os.remove(path)

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
y_settings = {'y0' : y0, 'y1' : y1, 'y2' : y2, 'ym' : ym, 'yf' : yf}
Din =   4.5     # lens entrance diameter
Dout =  2.4     # lens exit diameter
Dmax =  8.      # max diameter
D_settings = {'Din' : Din, 'Dout' : Dout, 'Dmax' : Dmax}
rIn =   0.006     # capillary radius
wall =   0.001 # |*50 make wider walls for structure visibility

# Hex structure parameters
nx_capillary = 13 
ny_bundle = 17

# Pinhole parameters
pinlen  = 0.01                # Length 
rpin    = 0.015               # Pinhole radius [mm] | =Dout/10. 
ypin    = 155 - pinlen        # Optical path position

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

def build_beamline(nrays=1e4):
    # Those parameters should be hel by some Lens object
    # FIXME - unfortunately they are used somewhere in screening ?
    beamLine = raycing.BeamLine(height=0)
    beamLine.y0 = y0
    beamLine.y1 = y1
    beamLine.y2 = y2
    beamLine.yf = yf
    beamLine.ym = ym
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

    # [1] - Create Lens object used for merging entrance
    # structure of capillaries and their shape in z direction 
    Lens = PolyCapillaryLens(y_settings = y_settings, \
                            D_settings = D_settings,\
                            material = mGlass)

    # Create object with (x,y) points describing hexagonal 
    # structure of polycapillary optics
    # This also should be a part of bigger Lens object
    entrance_Structure = HexStructure(nx_capillary=nx_capillary, \
                                    ny_bundle=ny_bundle, \
                                    rIn=rIn, wall=wall)

    # Self explanatory, TODO - pickle save|load this structure
    Lens.setStructure(entrance_Structure)

    # Show obtained structure and save as png
    entrance_Structure.test()

    # Total automatization of inserting capillaries
    beamLine = Lens.getCapillaries(beamLine)

    print 'Number of capillaries: ' + str(len(beamLine.capillaries))

    beamLine.exitScreen = rsc.Screen(beamLine,'ExitScreen', (0,y2,0))

    # Pinhole structure related container(s)
    beamLine.pinholes = []

    # Insert very short and very thin golden and lined capillaries 
    # into the focus, acting as a proper one way image sharpening pinhole
    # Focus size radius estimation
    focus_r = 2*rpin    # ? 

    for it in range(-10,11):
        x_in = it * focus_r

        print 'Inserting pinhole on z = 0, x =', str(x_in)
        # Pinholes
        pinhole = Pinhole(beamLine, 'pinh',\
#                roll = np.pi * np.cos(np.pi*it),\
                roll = 0,\
                x_in = x_in, r = rpin, y_in = ypin)
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
            rs.copy_beam(pinholeGlobalTotal, pinholeGlobal,\
                    good, includeState=True)

    # Expose screens to post pinholes beam
    postPinhole = scr.exposeScreens(beamLine, pinholeGlobalTotal,\
            [ypin, 200])

    # [2] - bypass the pinhole
    postNoPinhole = scr.exposeScreens(beamLine, beamCapillaryGlobalTotal,\
            [ypin, 200])

    # 
    outDict.update(prePinhole)
#    outDict.update(postPinhole)
    outDict.update(postNoPinhole)

    return outDict

rr.run_process = run_process

def main():
    beamLine = build_beamline()
    plot2D(beamLine)

    # Create xrtp.Plots in outside module  
    plots = scr.createPlots(beamLine, save=save)
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
    if save:
        plot.persistentName = 'pickle/' + plot.baseName + '.pickle'
    plots.append(plot)
    xrtr.run_ray_tracing(plots, repeats=repeats, beamLine=beamLine,\
            processes=processes)

    print 'Number of used capillaries: ', len(beamLine.capillaries)
#    return beamLine

if __name__ == '__main__':
    main()
