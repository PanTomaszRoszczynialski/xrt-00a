""" 
Copied from source viewer
hoping for easy modifications for testing various
object parameters and putting them directly between 
source and detector
"""

import numpy as np

import xrt.backends.raycing as raycing
import xrt.backends.raycing.sources as rs
#import xrt.backends.raycing.apertures as ra
import xrt.backends.raycing.oes as roe
import xrt.backends.raycing.run as rr
import xrt.backends.raycing.materials as rm
import xrt.plotter as xrtp
import xrt.runner as xrtr
import xrt.backends.raycing.screens as rsc

# Import custom shaped object
from customobjects import CustomShape

""" WELL DOCUMENTED PARAMETERS """
repeats = 1e4   # liczba
E0      = 9000  # [eV]
D_      = 40    # [mm] | constant screen position

xLimits = [-4.05, 4.05] # Plot limits
zLimits = xLimits       # axis square

processes = 4

""" GeometricSource():: PARAMETERS TO CHECK: """
bl_height   = 0.
bl_xzMax    = 0.

# x-direction parameters
distx       = 'flat'
dx          = 0.1
distxprime  = 'flat'
dxprime     = 0.1
# z-direction
distz       = 'flat'
dz          = 0.1
distzprime  = 'flat'
dzprime     = 0.1

def build_beamline(nrays=1000):
    beamLine = raycing.BeamLine(height=bl_height)

    # source appends itself to the provided beamline
    rs.GeometricSource(
        beamLine,'GeometricSource',(0,0,0), nrays=nrays,
        distx=distx, dx=dx, distxprime=distxprime, dxprime=dxprime,
        distz=distz, dz=dz, distzprime=distzprime, dzprime=dzprime,
        distE='lines', energies=(E0,), polarization='horizontal')

    # some mysterious parameters
    beamLine.xzMax = bl_xzMax

    # Container for multiple things
    beamLine.things = []

    # Set position and thickness
    limPhysY = [36.95, 39.95]
    limPhysX = [-3.0, 1.0]

    # Tested material
    mGold = rm.Material('Au', rho=19.3)
    mGlass  = rm.Material(('Si', 'O'), kind='mirror',\
                          quantities=(1, 2), rho=2.2)

#    thing = roe.OE(beamLine, 'thing', material=mGold,\
#                   shape=shape, limPhysY=limPhysY,\
#                   roll=0*np.pi/3., limPhysX=limPhysX)
    thing = CustomShape(beamLine, 'thing', material=mGold,\
                   R = 1.0, limPhysY=limPhysY,\
                   roll=0.*np.pi, limPhysX=limPhysX)

    # Contain
    beamLine.things.append(thing)

    # One screen is enough to see if object made any difference
    beamLine.myScreens = []
    beamLine.myScreens.append(rsc.Screen(\
                              beamLine,\
                              'd = {0:02d}'.format(D_),\
                              (0,D_,0)))
    
    return beamLine

def run_process(beamLine, shineOnly1stSource=False):
    beamSource = beamLine.sources[0].shine()

    # prepare empty python - dictionary for screens
    outDict = {}

    # Expose the thing
    # Global and local beams
    glob, loca = beamLine.things[0].multiple_reflect(beamSource)


    # Exposing screen to the beam
    outBeam = beamLine.myScreens[0].expose(glob)
    outDict['screen_{0:02d}'.format(D_)] = outBeam

    return outDict
rr.run_process = run_process

def main():
    beamLine = build_beamline()
    plots = []

    # Plot creation
    plot = xrtp.XYCPlot('screen_{0:02d}'.format(D_),(3,-1,),
        xaxis=xrtp.XYCAxis(r'$x$', 'mm',\
                           bins=256, ppb=2,\
                           limits=xLimits),
        yaxis=xrtp.XYCAxis(r'$z$', 'mm',\
                           bins=256, ppb=2,\
                           limits=zLimits),
        caxis='category', beamState='screen_{0:02d}'.format(D_),
        title='Distance from source = {0:02d} [mm]'.format(D_))

    # Names
    plot.baseName = 'object_over'
    plot.saveName = 'png/object/' + plot.baseName + '.png'

    plots.append(plot)

    # Plot creation
    plot = xrtp.XYCPlot('screen_{0:02d}'.format(D_),(1,),
        xaxis=xrtp.XYCAxis(r'$x$', 'mm',\
                           bins=256, ppb=2,\
                           limits=xLimits),
        yaxis=xrtp.XYCAxis(r'$z$', 'mm',\
                           bins=256, ppb=2,\
                           limits=zLimits),
        caxis=xrtp.XYCAxis('Reflections', 'number',\
                data=raycing.get_reflection_number,\
                bins=256, ppb=2, limits=None),
#        caxis='category', beamState='screen_{0:02d}'.format(D_),
        title='Distance from source = {0:02d} [mm]'.format(D_))

    # Names
    plot.baseName = 'object_out'
    plot.saveName = 'png/object/' + plot.baseName + '.png'

    plots.append(plot)

    xrtr.run_ray_tracing(plots, repeats=repeats,\
                         beamLine=beamLine,\
                         processes=processes)

if __name__ == '__main__':
    main()
