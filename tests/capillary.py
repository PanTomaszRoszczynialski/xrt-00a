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

class StraightCapillaryTest(object):
    """ Implement shining through a single pipe-like object """
    def __init__(self):
        """ Tania przestrzen reklamowa """

        # Default capillary position and radius
        self.x_in = 0.0
        self.r_in = 0.5

        # Default capillary dimensions:
        # 10cm capillary starting at 40mm
        self.y_start = 40
        self.y_end   = 140

        # Far screen distance from the end of capillary
        self.far_screen_dist = 20
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
                                     (0,self.y_end + self.far_screen_dist,0))

    def make_capillary(self):
        """ Not much here """
        cap_dims = [self.y_start, self.y_end]
        self.capillary = pl.StraightCapillary(self.beamLine,
                                              'straightcap',
                                              [0, 0, 0],
                                              x_in = self.x_in,
                                              r = self.r_in,
                                              roll = 0/3,
                                              limPhysY = cap_dims)

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

    def set_far_screen_distance(self, dist):
        """ """
        self.far_screen_dist = dist

    def set_capillary_position(self, pos):
        """ """
        self.x_in = pos

    def set_capillary_radius(self, rad):
        """ This is necessary """
        self.r_in = rad

    def set_capillary_length(self, val):
        """ Can't change the startpoint """
        if val is 0:
            print 'Capillary length too short'
            return
        self.y_end = self.y_start + val

    def set_prefix(self, fix):
        """ Przestrzen reklamowa """
        self.prefix = str(fix)
        print 'new prefix: ', fix

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

    def run_it(self):
        """ Prepares plots """

        self.make_it()

        # xrt.Plot constants
        bins = 256
        plots = []

        # Screen at the end of capillary
        limits_near_x = [self.x_in - 1.2 * self.r_in,
                         self.x_in + 1.2 * self.r_in]
        limits_near_z = [-1.2 * self.r_in,
                         1.2 * self.r_in]
        plot = xrtp.XYCPlot(
            # Using named parameters might be good here TODO
            'ExitScreen', (1, 3,),
            beamState='ExitScreen',
            # This is still inside plot
            xaxis=xrtp.XYCAxis(r'$x$', 'mm',
                               bins=bins,
                               ppb=2,
                               limits=limits_near_x),
            # As is this
            yaxis=xrtp.XYCAxis(r'$z$', 'mm',
                               bins=bins,
                               ppb=2,
                               limits=limits_near_z),
            # Colorbar and sidebar histogram
            caxis=xrtp.XYCAxis('Reflections',
                               'number',
                               data=raycing.get_reflection_number,
                               bins=bins,
                               ppb=2,
                               limits=[0,10])
            ) # plot ends here
        plot.title = 'Capillary position'
        plot.saveName = 'png/tests/near/' + self.prefix + '_near.png'
        plots.append(plot)

        # Plot far after the exit
        limits_far = [-5, 5]
        plot = xrtp.XYCPlot(
            'FarScreen', (1, 3,),
            xaxis=xrtp.XYCAxis(r'$x$', 'mm',
                               bins=bins,
                               ppb=2,
                               limits=limits_far),
            yaxis=xrtp.XYCAxis(r'$z$', 'mm',
                               bins=bins,
                               ppb=2,
                               limits=limits_far),
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
        plot.saveName = 'png/tests/far/' + self.prefix + '_far.png'
        plots.append(plot)

        xrtr.run_ray_tracing(plots, repeats=20, beamLine=self.beamLine,\
                processes=1)

def run_test():
    """ Full test """
    test = StraightCapillaryTest()

    # Define values to be tested
    radiuses  = [0.1, 0.2, 0.5, 0.7, 1.0]
    positions = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    lengths   = [10, 30, 50, 80, 100, 120, 150, 200]
    radiuses = [0.5]
    positions = [0.0]
    lengths = [10]
    screen_dsits = [60, 70, 80, 90, 100]

    # Run
    for radius in radiuses:
        test.set_capillary_radius(radius)
        for position in positions:
            test.set_capillary_position(position)
            for length in lengths:
                test.set_capillary_length(length)
                for dist in screen_dsits:
                    test.set_far_screen_distance(dist)
                    fix = 'rad_' + str(radius)
                    fix += '___pos_' + str(position)
                    fix += '___len_' + str(1000+length)
                    fix += '___fsd_' + str(1000+dist)
                    test.set_prefix(fix)
                    test.run_it()

class TaperedCapillaryTest(object):
    """ This class is supposed to help with testing more complex
    capillary shapes: with straight axis and varying radius """
    def __init__(self):
        """ Tania przestrzen reklamowa """

        # Set start and end of a capillary
        self.y_start = 40
        self.y_end = 140
        self.x_in = 0.0

        # Set constant radius
        self.r_in = 0.5
        self.r_out = 0.5
