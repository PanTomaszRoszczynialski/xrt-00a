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
    
mGlass = rm.Material(('Si', 'O'), quantities=(1, 2), rho=2.2)

repeats = 6*1500 # number of ray traycing iterations
E0 = 9000.
rSample = 120 # starting position of the lens
f = rSample + 150 # y length in mm from foucs to the end of the lens
screen1_pos = rSample + 100 
screen2_pos = f + 5 # distance @vincze == 10cm
max_plots = 0
r0 = 0.002
rOut = 0.002
wall = 0.0005
plot2D_yLim = [-0.05, 0.05]
plot_main_lim = 0.45 # min 2*r0 for capillary entrance imaging
layers = 10 # number of hexagonal layers
nRefl = 11 # number of reflections
nReflDisp = 12 # unused
xzPrimeMax = 3.
# Pickle saving: None for no saving
persistentName=None #'phase_space__energy.pickle'
# some fun parameter
k_ = 256.0

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
        return -k_*np.cosh((s-75.)/75./k_) +  k_*np.cosh(1/k_)
        

    def local_x0Prime(self, s):
        return -1.0/75.0 * np.sinh((s-75.0)/75.0/k_)

    def local_r0(self, s):  # radius of capillary (s)
#        return self.ar * (s-self.s0)**2 + self.br
        return -self.ar *(s-self.s0) + self.br

    def local_r0Prime(self, s):
#        return self.ar * 2 * (s-self.s0)
        return -self.ar

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
        beamLine, 'GeometricSource', (0,0,0), nrays=nrays,
        dx=0., dz=0., distxprime='annulus',
        distE='lines', energies=(E0,), polarization='horizontal')                 
    # yo    
    beamLine.fsm1 = rsc.Screen(beamLine, 'DiamondFSM1', (0,screen1_pos,0))
    
    # try to remove superflous container
    #beamLine.capillaries = []
    beamLine.xzMax = 0 # no ide what this does
    # this parameter should be @line 8
    alpha = 0.000 # hopefully milliradian
    roll = 0 # test if this rotates whole object
    capillary = StraightCapillary(
        beamLine, 'StraightCapillary', [0,0,0], roll=roll,
        material=mGlass, limPhysY=[rSample*np.cos(alpha), f],
        order=8, f=f, rSample=rSample, entranceAlpha=alpha, rIn=r0, rOut=rOut)
    beamLine.capillary = capillary         
#    beamLine.capillaries.append(capillary)         
    
    if beamLine.xzMax < capillary.b0:
        beamLine.xzMax = capillary.b0
    beamLine.xzMax += 2*r0
    
    n=2    # one layer..
    beamLine.sources[0].dxprime = 0, np.arcsin((2*n+1) * (r0+wall) / rSample)
    beamLine.fsm2 = rsc.Screen(beamLine,'DiamondFSM2', (0,screen2_pos,0))
    beamLine.myFsms = []
    for it in range(0,max_plots):
        beamLine.myFsms.append(rsc.Screen(beamLine,'myScreen{0:02d}'.format(it),(0,rSample+145*it,0)))

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

    xLimits = [-r0*1.4, 1.4*r0]
    xpLimits = [-0.2, 0.2]
    zLimits = [-r0*1.4, 1.4*r0]
#    yLimits=None
    cLimits = None #[8900,9100]
    # at the entrance
    """
    PHASE SPACE PLOT
    """
    plot = xrtp.XYCPlot('beamFSM2', (1,),
        xaxis=xrtp.XYCAxis(r"$z$", 'mm', data=raycing.get_z, bins=256, ppb=2, limits=zLimits),
        yaxis=xrtp.XYCAxis(r"$z'$", 'mrad', data=raycing.get_zprime, bins=256, ppb=2, limits=xpLimits),
#        caxis='category', 
        caxis=xrtp.XYCAxis("Phse shift", 'mrad',data=raycing.get_phase_shift, bins=256, ppb=2, limits=cLimits),
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
    plot = xrtp.XYCPlot('beamFSM2', (1,),
        xaxis=xrtp.XYCAxis(r"$x$", 'mm', data=raycing.get_x, bins=256, ppb=2, limits=xLimits),
        yaxis=xrtp.XYCAxis(r"$z$", 'mm', data=raycing.get_z, bins=256, ppb=2, limits=zLimits),
#        caxis='category', 
        caxis=xrtp.XYCAxis("Reflections", 'nr of',data=raycing.get_reflection_number, bins=256, ppb=2, limits=[0,10]),
        beamState='beamFSM2', title='Real Space', aspect='auto',
        persistentName=persistentName)
    # setting persistentName saves data into a python pickle, and might be
    # unhealthy if pickle isn't cleared/deleted when plotted data changes
    plot.baseName = 'realSpace'
    plot.saveName = plot.baseName + '.png'    
    plots.append(plot)
    
    # ITERATING OVER PLOTS {}
    xrtr.run_ray_tracing(plots, repeats=repeats, beamLine=beamLine, processes=1)
    
    # savemat() takes a dict of names later loaded into matlab and objects
    # we want to save,
    scipy.io.savemat('Phase_xrt_data.mat',{'Phase_total2D_RGB':plots[0].total2D_RGB,
                                 'Phase_total2D':plots[0].total2D,
                                 'Phase_caxis_total':plots[0].caxis.total1D,
                                 'Phase_caxis_total_RGB':plots[0].caxis.total1D_RGB})
    scipy.io.savemat('RSpace_xrt_data.mat',{'RSpace_total2D_RGB':plots[1].total2D_RGB,
                                 'RSpace_total2D':plots[1].total2D,
                                 'RSpace_caxis_total':plots[1].caxis.total1D,
                                 'RSpace_caxis_total_RGB':plots[1].caxis.total1D_RGB})                                 
    # just for debug 
    return plot                                 
    
    
if __name__ == '__main__':
#    PlotMono.plot2D(build_beamline(),f)
    main()