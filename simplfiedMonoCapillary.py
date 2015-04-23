# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 10:51:55 2015

@author: Vladimir Putin
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import PlotMono

import xrt.backends.raycing as raycing
import xrt.backends.raycing.sources as rs
#import xrt.backends.raycing.apertures as ra
import xrt.backends.raycing.oes as roe
import xrt.backends.raycing.run as rr
import xrt.backends.raycing.materials as rm
import xrt.plotter as xrtp
import xrt.runner as xrtr
import xrt.backends.raycing.screens as rsc

# see XYCAxis constructor:
#from xrt.backends import raycing

# for saving into .mat file
import scipy.io
    
# ray traycing settings    
mGlass = rm.Material(('Si', 'O'), quantities=(1, 2), rho=2.2)
repeats = 1e4    # number of ray traycing iterations
E0 = 9000.          # energy in electronoVolts
nRefl = 30          # number of reflections

# capillary shape parameters
rSample = 30.0 # starting position of the lens
L_ = 100.0
f = rSample + L_ # y length in mm from foucs to the end of the lens
r0 = 0.02
rOut = 0.02
wall = 0.005

# parameters for local_x0 function for actual shape definition
y_in    = 1.1            # entrance height
rS      = float(rSample)    # light source - capillary distance 
# Cosh parameter for tangential ray entrance
a_      = -L_/2/np.arcsinh(-y_in/rS)
print a_, y_in/rS

# image acquisition
screen1_pos = rSample     # not really used
screen2_pos = f + 1             # first image position outside capillary
max_plots = 0                   # for imaging different position at once| 0=off

# Pickle saving: None for no saving
persistentName=None #'phase_space__energy.pickle'

# Source parameters
distx       = 'flat'
dx          = 0.1
distxprime  = 'normal'
dxprime     = 0.1
# z-direction
distz       = 'flat'
dz          = 0.1
distzprime  = 'normal'
dzprime     = 0.1


class StraightCapillary(roe.OE):
    def __init__(self, *args, **kwargs):
        self.rSample = kwargs.pop('rSample')
        self.entranceAlpha = kwargs.pop('entranceAlpha')
        self.f = kwargs.pop('f') # 
        self.r0in = kwargs.pop('rIn')
        self.r0out = kwargs.pop('rOut')
        roe.OE.__init__(self, *args, **kwargs)

        s0 = self.f - self.rSample * np.cos(self.entranceAlpha)
        self.a0 = -np.tan(self.entranceAlpha) / 2 / s0
        self.b0 = 0.5*self.rSample * np.sin(self.entranceAlpha) - self.a0 * s0**2
        self.b0 = 0.
        self.s0 = s0
        self.ar = (self.r0out-self.r0in) / s0
        self.br = self.r0in
        self.isParametric = True

    def local_x0(self, s):  # axis of capillary, x(s)
        # s*0 is needed for this method to act as a function rather than variable?
        return -a_*np.cosh((s-L_/2.)/a_) +  a_*np.cosh(L_/2./a_) + y_in
        

    def local_x0Prime(self, s):
        return -np.sinh((s-L_/2.)/a_)

    def local_r0(self, s):  # radius of capillary (s)
#        return self.ar * (s-self.s0)**2 + self.br
        return -self.ar *(s-self.s0) + self.br #self.br*(2+np.cos(np.pi/2*(s-L_/2.)/L_/2.))

    def local_r0Prime(self, s):
#        return self.ar * 2 * (s-self.s0)
        return -self.ar #+ self.br * ( - np.sin(np.pi/2*(s-L_/2.)/L_/2.))*np.pi/2/L_/2.

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
        
def build_beamline(nrays=500):
    beamLine = raycing.BeamLine(height=0)
    # source checked @SourceView.py
    rs.GeometricSource(
        beamLine, 'GeometricSource', (0,0,0), nrays=nrays,
        distx=distx, dx=dx, distxprime=distxprime, dxprime=dxprime,
        distz=distz, dz=dz, distzprime=distzprime, dzprime=dzprime,
        distE='lines', energies=(E0,), polarization='horizontal')            
    
#    rs.GeometricSource(
#        beamLine, 'GeometricSource', (0,0,0), nrays=nrays,
#        dx=0., dz=0., distxprime='annulus',
#        distE='lines', energies=(E0,), polarization='horizontal')                
#    
    # yo    
    beamLine.fsm1 = rsc.Screen(beamLine, 'DiamondFSM1', (0,screen1_pos,0))
    
    alpha = np.arctan(y_in/rSample)
    alpha = 0
    roll = 0 # test if this rotates whole object
    limPhysY=[0, f]
#    limPhysY=[0, f]

    
    capillary = StraightCapillary(
        beamLine, 'StraightCapillary', [0,0,0], roll=roll,
        material=mGlass, limPhysY=limPhysY,
        order=8, f=f, rSample=rSample, entranceAlpha=alpha, rIn=r0, rOut=rOut)
    beamLine.capillary = capillary         
#    beamLine.capillaries.append(capillary)         
    
    beamLine.fsm2 = rsc.Screen(beamLine,'DiamondFSM2', (0,screen2_pos,0))
    beamLine.myFsms = []
    for it in range(0,max_plots):
        beamLine.myFsms.append(rsc.Screen(beamLine,
                                          'myScreen{0:02d}'.format(it),
                                          (0,f +rSample + rSample/10*it,0)))

    return beamLine
         
def run_process(beamLine, shineOnly1stSource=False):
    beamSource = beamLine.sources[0].shine()
    # at the entrance
    beamFSM1 = beamLine.fsm1.expose(beamSource)
    outDict = {'beamSource': beamSource, 'beamFSM1': beamFSM1}
    beamCapillaryGlobalTotal = None
    capillary = beamLine.capillary
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
    outDict['myBeam_after_local'] = beamCapillaryLocalN
    outDict['myBeam_after_global'] = beamCapillaryGlobalTotal

    # Create second screen
    beamFSM2 = beamLine.fsm2.expose(beamCapillaryGlobalTotal)
    outDict['beamFSM2'] = beamFSM2
    beamFsms = []
    for it in range(0,max_plots):
        beamFsms.append(beamLine.myFsms[it].expose(beamCapillaryGlobalTotal))
        outDict['myExposedScreen{0:02d}'.format(it)] = beamFsms[it]
    return outDict
    
rr.run_process = run_process  


def main():
    beamLine = build_beamline()
    plots = []

#    xLimits = [y_in - 3.6*r0, y_in + 3.6*r0]
    xpLimits = [-0.55, 0.55]
#    zLimits = xLimits
    zLimits = [-r0*3.6, 3.6*r0]
    xLimits = None
    zLimits = xLimits
#    yLimits=None
    cLimits = [-3,3] #[8900,9100]
    # at the entrance
    """
    PHASE SPACE PLOT
    """
    plot = xrtp.XYCPlot('beamFSM2', (1,3),
        xaxis=xrtp.XYCAxis(r"$z$", 'mm', data=raycing.get_z, bins=256, ppb=2, limits=zLimits),
        yaxis=xrtp.XYCAxis(r"$z'$", 'mrad', data=raycing.get_zprime, bins=256, ppb=2, limits=xpLimits),
#        caxis='category', 
        caxis=xrtp.XYCAxis("Phse shift", 'rad',data=raycing.get_phase_shift, bins=256, ppb=2, limits=cLimits),
        beamState='beamFSM2', title='Phase Space', aspect='auto',
        persistentName=persistentName)
    # setting persistentName saves data into a python pickle, and might be
    # unhealthy if pickle isn't cleared/deleted when plotted data changes
    plot.baseName = 'phaseSpace'
    plot.saveName = plot.baseName + '.png'
    plots.append(plot)
    
    """
    REAL SPACE PLOT
    """
    plot = xrtp.XYCPlot('beamFSM2', (1,3),
        xaxis=xrtp.XYCAxis(r"$x$", 'mm', data=raycing.get_x, bins=256, ppb=2, limits=xLimits),
        yaxis=xrtp.XYCAxis(r"$z$", 'mm', data=raycing.get_z, bins=256, ppb=2, limits=zLimits),
#        caxis='category', 
        caxis=xrtp.XYCAxis("Reflections", 'nr of',data=raycing.get_reflection_number, bins=256, ppb=2, limits=[0,nRefl]),
        beamState='beamFSM2', title='Real Space', aspect='auto',
        persistentName=persistentName)
    # setting persistentName saves data into a python pickle, and might be
    # unhealthy if pickle isn't cleared/deleted when plotted data changes
    plot.baseName = 'realSpace'
    plot.saveName = plot.baseName + '.png'    
    plots.append(plot)
    
    # ITERATING OVER PLOTS {}
    for it in range(0,max_plots):
        plot = xrtp.XYCPlot('myExposedScreen{0:02d}'.format(it), (1,),
            xaxis=xrtp.XYCAxis(r'$x$', 'mm', bins=256, ppb=2, limits=xLimits),
            yaxis=xrtp.XYCAxis(r'$z$', 'mm', bins=256, ppb=2, limits=zLimits),
            caxis='category', beamState='myExposedScreen{0:02d}'.format(it), title=str(it))
        plot.baseName = 'thin_cap_dist_' + str(110+it)
        plot.saveName = plot.baseName + '.png'
        plots.append(plot)    
    xrtr.run_ray_tracing(plots, repeats=repeats, beamLine=beamLine, processes=7)
    
    # savemat() takes a dict of names later loaded into matlab and objects
    # we want to save,
#    scipy.io.savemat('Phase_xrt_data.mat',{'Phase_total2D_RGB':plots[0].total2D_RGB,
#                                 'Phase_total2D':plots[0].total2D,
#                                 'Phase_caxis_total':plots[0].caxis.total1D,
#                                 'Phase_caxis_total_RGB':plots[0].caxis.total1D_RGB})
#    scipy.io.savemat('RSpace_xrt_data.mat',{'RSpace_total2D_RGB':plots[1].total2D_RGB,
#                                 'RSpace_total2D':plots[1].total2D,
#                                 'RSpace_caxis_total':plots[1].caxis.total1D,
#                                 'RSpace_caxis_total_RGB':plots[1].caxis.total1D_RGB})                                 
    # just for debug 
    return plot                                 
    
    
if __name__ == '__main__':
#    PlotMono.plot2D(build_beamline(),f)
    main()