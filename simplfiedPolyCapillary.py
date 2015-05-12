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

import xrt.backends.raycing.apertures as ra


# for saving into .mat file
import scipy.io
    
# ray traycing settings    
mGlass = rm.Material(('Si', 'O'), quantities=(1, 2), rho=2.2)
repeats = 1e4       # number of ray traycing iterations
E0 = 9000.          # energy in electronoVolts
nRefl = 30         # number of reflections

# capillary shape parameters
rSample = 30.0              # starting position of the lens
L_      = 200.0               # length of the lens
f       = rSample + L_     # y length in mm from foucs to the end of the lens
r0 = 0.002*10
rOut = 0.002*10
wall = 0.0005

# parameters for local_x0 function for actual shape definition
y_in    = 0.12             # entrance height
rS      = float(rSample)    # light source - capillary distance 
# Cosh parameter for tangential ray entrance
a_      = -L_/2.0/np.arcsinh(-y_in/rS)
print a_, y_in/rS

# image acquisition
screen1_pos = rSample     # not really used
screen2_pos = f + 0             # first image position outside capillary
max_plots = 0                   # for imaging different position at once| 0=off

# Pickle saving: None for no saving
persistentName = 'pickle/polyCapExit.pickle' #'realSpae.pickle' 

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
        return -a_*np.cosh((s-L_/2.0)/a_) +  a_*np.cosh(L_/2.0/a_) + y_in
        

    def local_x0Prime(self, s):
        return -np.sinh((s-L_/2.0)/a_)

    def local_r0(self, s):  # radius of capillary (s)
        return self.ar *(s-self.s0)**2 + self.br #+ self.br*(2.0+np.cos(np.pi/2.0*(s-L_/2.0)/L_/2.0))

    def local_r0Prime(self, s):
        return self.ar * 2 * (s-self.s0) #+ self.br * ( - np.sin(np.pi/2.0*(s-L_/2.0)/L_/2.0))*np.pi/2.0/L_/2.0


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
        
def build_beamline(nrays=1e4):
    beamLine = raycing.BeamLine(height=0)
#    rs.GeometricSource(
#        beamLine, 'GeometricSource', (0,0,0), nrays=nrays,
#        dx=0., dz=0., distxprime='annulus',
#        distE='lines', energies=(E0,), polarization='horizontal')    
    rs.GeometricSource(
        beamLine,'GeometricSource',(0,0,0), nrays=nrays,
        distx=distx, dx=dx, distxprime=distxprime, dxprime=dxprime,
        distz=distz, dz=dz, distzprime=distzprime, dzprime=dzprime,
        distE='lines', energies=(E0,), polarization='horizontal')         
    # yo    
    beamLine.fsm1 = rsc.Screen(beamLine, 'DiamondFSM1', (0,screen1_pos,0))
    
    # try to remove superflous container
    beamLine.capillaries = []

    alpha = 0.000   # hopefully milliradian
    N_ = 12
    for it in range(N_):
        roll = it*2*np.pi/N_
        capillary = BentCapillary(
            beamLine, 'BentCapillary', [0,0,0], roll=roll,
            material=mGlass, limPhysY=[rSample*np.cos(alpha), f],
            order=8, f=f, rSample=rSample, entranceAlpha=alpha, rIn=r0, rOut=rOut)
        # 
        beamLine.capillaries.append(capillary)         
    
    # debug cout
    print "Number of capillaries: " + str(len(beamLine.capillaries))
    
    # prepare screen fo all of them
    beamLine.fsm2 = rsc.Screen(beamLine,'DiamondFSM2', (0,screen2_pos,0))
    beamLine.myFsms = []
    
    # Iterate from exit to focus (symmetric atm), save distances for names
    beamLine.myScreens_pos = []
    for it in range(1,max_plots+1):
        tmp_pos = f + it*rSample/max_plots
        beamLine.myScreens_pos.append(tmp_pos)
        beamLine.myFsms.append(rsc.Screen(beamLine,
                                          'myScreen{0:02d}'.format(it),
                                          (0, tmp_pos, 0)))
                                          
    # Insert slit                                          
    slitDx = 0.2
    slitDz = 0.1                                     
    beamLine.slit = ra.RectangularAperture(
        beamLine, 'squareSlit', [0, f + 2.2*rSample ,0], ('left', 'right', 'bottom', 'top'),
        [-slitDx, slitDx, -slitDz, slitDz])
        
    # Insert screen after slit
    beamLine.sltScreen = rsc.Screen(beamLine,'DiamondFSM2', (0,f+2.5*rSample,0))        

    return beamLine
         
def run_process(beamLine, shineOnly1stSource=False):

    beamSource = beamLine.sources[0].shine()
    # at the entrance | unused
    beamFSM1 = beamLine.fsm1.expose(beamSource)
    outDict = {'beamSource': beamSource, 'beamFSM1': beamFSM1}
    
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
    beamFSM2 = beamLine.fsm2.expose(beamCapillaryGlobalTotal)
    outDict['beamFSM2'] = beamFSM2
    
    # For future use | this for should probably end before the slit
    # and another one should be run after | ?
    beamFsms = []
    for it in range(0,max_plots): # 0
        beamFsms.append(beamLine.myFsms[it].expose(beamCapillaryGlobalTotal))
        outDict['myExposedScreen{0:02d}'.format(it)] = beamFsms[it]

    # Push photons through the slit
    beamPastSlit, dupa = beamLine.slit.propagate(beamCapillaryGlobalTotal,needNewGlobal=True)
    
    # Observe what's behind
    pastSlitScreen = beamLine.sltScreen.expose(beamPastSlit)
    outDict['slitScreen'] = pastSlitScreen
    return outDict
    
rr.run_process = run_process  


def main():
    beamLine = build_beamline()
    plots = []

    limit_r = y_in + 1.6 * r0     
    xLimits = [- limit_r, limit_r]
    zLimits = [-limit_r, limit_r]

    
    """
    REAL SPACE PLOT
    """
    plot = xrtp.XYCPlot('beamFSM2', (1,),
        xaxis=xrtp.XYCAxis(r"$x$", 'mm', data=raycing.get_x, bins=256, ppb=2, limits=xLimits),
        yaxis=xrtp.XYCAxis(r"$z$", 'mm', data=raycing.get_z, bins=256, ppb=2, limits=zLimits),
#        caxis='category', 
        caxis=xrtp.XYCAxis("Reflections", 'num of',data=raycing.get_reflection_number, bins=256, ppb=2, limits=[0, nRefl]),
        beamState='beamFSM2', title='Detector at exit', aspect='auto',
        persistentName=persistentName)
    # setting persistentName saves data into a python pickle, and might be
    # unhealthy if pickle isn't cleared/deleted when plotted data changes
    plot.baseName = 'Detector_at__' + str(f)
    plot.saveName = 'png/' + plot.baseName + '.png'    
    plots.append(plot)
    
    """
    BEHIND THE SLIT
    """
    limit_r = 0.5
    xLimits = [- limit_r, limit_r]
    zLimits = [-limit_r, limit_r]
    plot = xrtp.XYCPlot('slitScreen', (1,),
        xaxis=xrtp.XYCAxis(r"$x$", 'mm', data=raycing.get_x, bins=256, ppb=2, limits=xLimits),
        yaxis=xrtp.XYCAxis(r"$z$", 'mm', data=raycing.get_z, bins=256, ppb=2, limits=zLimits),
        caxis='category', 
#        caxis=xrtp.XYCAxis("Reflections", 'num of',data=raycing.get_reflection_number, bins=256, ppb=2, limits=[0, nRefl]),
        beamState='slitScreen', title='After the slit', aspect='auto',
        persistentName=None)
    plot.baseName = 'slit_0A'
    plot.saveName = 'png/' + plot.baseName + '.png'    
    plots.append(plot)
    
    # ITERATING OVER PLOTS {}
    cLimits = [0, nRefl]
    for it in range(0,max_plots):
        tmp_name = 'Detector_at_' + str(beamLine.myScreens_pos[it])
        plot = xrtp.XYCPlot('myExposedScreen{0:02d}'.format(it), (1,3),
            xaxis=xrtp.XYCAxis(r'$x$', 'mm', bins=256, ppb=2, limits=xLimits),
            yaxis=xrtp.XYCAxis(r'$z$', 'mm', bins=256, ppb=2, limits=zLimits),
            caxis=xrtp.XYCAxis('Reflections', 'number', data=raycing.get_reflection_number, bins=256, ppb=2, limits=cLimits),
            beamState='myExposedScreen{0:02d}'.format(it), title=tmp_name)
        plot.baseName = tmp_name
        plot.saveName = 'png/' + plot.baseName + '.png'
        plot.persistentName = 'pickle/' + plot.baseName + '.pickle'
        plots.append(plot)    
    xrtr.run_ray_tracing(plots, repeats=repeats, beamLine=beamLine, processes=1)
    
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
    # just for debug x`
    return plot                                 
    
if __name__ == '__main__':
#    PlotMono.plot2D(build_beamline(),f)
    main()