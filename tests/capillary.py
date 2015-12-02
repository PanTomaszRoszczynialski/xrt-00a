import matplotlib as mpl
mpl.use('Agg')

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

        self.x_in = x_in
        self.r_in = r_in

        # Default capillary dimensions:
        # 10cm capillary starting at 40mm
        self.y_start = 40
        self.y_end   = 140

        # Savename prefix (colon helps with vim-folding)
        self.prefix = 'dupa';

    def make_it(self):
        """ Run this everytime any of the parameters are re-set """
        self.make_source()
        self.make_capillary()
        self.make_screens()
        self.make_run_process()

    def make_screens(self):
        """ One screen at the exit and one somewhere far """
        self.beamLine.exitScreen = rsc.Screen(self.beamLine,
                                             'Exitscreen',
                                             (0, self.y_end, 0))

        self.beamLine.farScreen = rsc.Screen(self.beamLine,
                                             'FarScreen',
                                             (0,self.y_end+100,0))

    def make_capillary(self):
        """ Not much here """
        self.capillary = pl.StraightCapillary(self.beamLine,
                                              'straightcap',
                                              [0, 0, 0],
                                              x_in = self.x_in,
                                              r = self.r_in,
                                              roll = 0/3,
                                              limPhysY = [40, 142])

    def make_source(self):
        """ Many source parameters are set here by hand """

        # Source parameters
        nrays       = 2000
        E0          = 9000
        # x-direction
        distx       = 'flat'
        dx          = self.r_in * 1.5
        distxprime  = 'flat'
        dxprime     = 0.02
        # z-direction
        distz       = 'flat'
        dz          = self.r_in * 1.5
        distzprime  = 'flat'
        dzprime     = 0.02

        # Sources appends themselves to the beamline
        # so we recreate the beamline to always use bl.sources[0]
        self.beamLine = raycing.BeamLine(height=0)
        ss.FittedSource(
            self.beamLine,'DirectedSource',(0,39.9,0), nrays=nrays,
            distx=distx, dx=dx, distxprime=distxprime, dxprime=dxprime,
            distz=distz, dz=dz, distzprime=distzprime, dzprime=dzprime,
            distE='lines', energies=(E0,), polarization='horizontal')

    def set_capillary_length(self, val):
        """ Can't change the startpoint """
        self.y_end = self.y_start + val

    def set_prefix(self, fix):
        """ Przestrzen reklamowa """
        self.prefix = str(fix)

    def make_run_process(self):
        """ Overloads xrt method for photon generation """
        def local_process(beamLine, shineOnly1stSource=False):
            beamTotal = None
            x, z = self.capillary.entrance_point()
            hitpoint = (x, 40, z)
            beamSource = beamLine.sources[0].shine(hitpoint=hitpoint)

            beamTotal, _ = self.capillary.multiple_reflect(beamSource,\
                                              maxReflections=8)

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

        self.make_it()

        beamLine = self.beamLine

        # xrt.Plot constants
        bins = 256
        plots = []

        # Screen at the end of capillary
        limits_near = [-1.2 * self.r_in, 1.2 * self.r_in]
        plot = xrtp.XYCPlot(
            # Using named parameters might be good here TODO
            'ExitScreen', (1, 3,),
            beamState='ExitScreen',
            # This is still inside plot
            xaxis=xrtp.XYCAxis(r'$x$', 'mm',
                               bins=bins,
                               ppb=2,
                               limits=limits_near),
            # As is this
            yaxis=xrtp.XYCAxis(r'$z$', 'mm',
                               bins=bins,
                               ppb=2,
                               limits=limits_near),
            # Colorbar and sidebar histogram
            caxis=xrtp.XYCAxis('Reflections',
                               'number',
                               data=raycing.get_reflection_number,
                               bins=bins,
                               ppb=2,
                               limits=[0,10])
            ) # plot ends here
        plot.title = 'Capillary position'
        plot.saveName = 'png/tests/' + self.prefix + '_near.png'
        plots.append(plot)

        # Plot far after the exit
        plot = xrtp.XYCPlot(
            'FarScreen', (1, 3,),
            xaxis=xrtp.XYCAxis(r'$x$', 'mm',
                               bins=bins,
                               ppb=2,
                               limits=None),
            yaxis=xrtp.XYCAxis(r'$z$', 'mm',
                               bins=bins,
                               ppb=2,
                               limits=None),
            beamState='FarScreen',
            # Colorbar and sidebar histogram 
            caxis=xrtp.XYCAxis('Reflections',
                               'number',
                               data=raycing.get_reflection_number,
                               bins=bins,
                               ppb=2,
                               limits=[0,10])
            ) # different plot ends here
        plot.title = 'Divergence observations'
        plot.saveName = 'png/tests/' + self.prefix + '_far.png'
        plots.append(plot)

        xrtr.run_ray_tracing(plots, repeats=20, beamLine=beamLine,\
                processes=1)

def run_test():
    """ Full test """
    x_in = 0.1
    r_in = 0.5
    test = CapillaryTest(x_in, r_in)
    for it in range(20):
        newlen = 10 + 20 * it
        print 'Current length:', newlen, '[mm]'
        prefix = 'capillary_length__' + str(1000+newlen)
        test.set_prefix(prefix)
        test.set_capillary_length(newlen)
        test.run_test()
