import numpy as np
import xrt.runner as xrtr
import xrt.backends.raycing.screens as rsc
import xrt.plotter as xrtp
import xrt.backends.raycing.run as rr
import xrt.backends.raycing as raycing

from lenses import polycapillary as pl
from sources import special_sources as ss

# Create source
# Shine on a capillary
# show on screen
# done!

class CapillaryTest(object):
    """ Implement shining through a single pipe-like object """
    def __init__(self, x_in, r_in):
        """ Take entrance point and radius of the capillary """
        # xrt commons
        self.beamLine = raycing.BeamLine(height=0)

        # Source parameters
        nrays       = 300
        E0          = 9000
        # x-direction
        distx       = 'flat'
        dx          = r_in * 1.5
        distxprime  = 'flat'
        dxprime     = 0.02
        # z-direction
        distz       = 'flat'
        dz          = r_in * 1.5
        distzprime  = 'flat'
        dzprime     = 0.02

        ss.FittedSource(
            self.beamLine,'DirectedSource',(0,40,0), nrays=nrays,
            distx=distx, dx=dx, distxprime=distxprime, dxprime=dxprime,
            distz=distz, dz=dz, distzprime=distzprime, dzprime=dzprime,
            distE='lines', energies=(E0,), polarization='horizontal')

        args = [self.beamLine, 'bent', [0,0,0]]
        self.capillary = pl.StraightCapillary(self.beamLine,
                                              'straight',
                                              [0, 0, 0],
                                              x_in = 0.5,
                                              r = 1.1,
                                              roll = 1/3,
                                              limPhysY = [40, 140])

        # One screen at the exit and one somewhere far
        self.beamLine.exitScreen = rsc.Screen(self.beamLine,
                                             'EntranceScreen',
                                             (0, 140, 0))

        self.beamLine.farScreen = rsc.Screen(self.beamLine,
                                             'FarScreen',
                                             (0,200,0))
    def make_run_process(self):
        """ Overloads xrt method for photon generation """
        def local_process(beamLine, shineOnly1stSource=False):
            beamTotal = None
            x, z = self.capillary.entrance_point()
            hitpoint = (x, 40, z)
            beamSource = beamLine.sources[0].shine(hitpoint=hitpoint)

            beamTotal, dupa = self.capillary.multiple_reflect(beamSource,\
                                              maxReflections=10)

            # Expose screens
            exitScreen = beamLine.exitScreen.expose(beamTotal)
            farScreen = beamLine.farScreen.expose(beamTotal)

            # Show beamlines after exposition
            outDict = {'ExitScreen' : exitScreen}
            outDict['FarScreen'] = farScreen

            return outDict

        rr.run_process = local_process

    def run_test(self):
        """ Prepares plots """

        beamLine = self.beamLine

        # xrt.Plot constants
        bins = 256
        # limits = [-20,20]
        limits = None
        plots = []

        # Screen at the end of capillary
        plot = xrtp.XYCPlot(
            # Using named parameters might be good here TODO
            'ExitScreen', (1, 3,),
            caxis='category', beamState='ExitScreen',
            # This is still inside plot
            xaxis=xrtp.XYCAxis(r'$x$', 'mm',
                               bins=bins,
                               ppb=2,
                               limits=limits),
            # As is this
            yaxis=xrtp.XYCAxis(r'$z$', 'mm',
                               bins=bins,
                               ppb=2,
                               limits=limits)
            ) # plot ends here
        plot.title = 'Capillary position'
        plots.append(plot)

        plot = xrtp.XYCPlot(
            'FarScreen', (1, 3,),
            xaxis=xrtp.XYCAxis(r'$x$', 'mm',
            bins=bins, ppb=2, limits=None),
            yaxis=xrtp.XYCAxis(r'$z$', 'mm',
            bins=bins, ppb=2, limits=None),
            caxis='category', beamState='FarScreen')
        plot.title = 'Divergence observations'
        plots.append(plot)

        xrtr.run_ray_tracing(plots, repeats=100, beamLine=beamLine,\
                processes=1)

def run_test():
    """ Full test """
    x_in = 2
    r_in = 0.3
    test = CapillaryTest(x_in, r_in)
    test.make_run_process()
    test.run_test()

