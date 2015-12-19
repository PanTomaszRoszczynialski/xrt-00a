import matplotlib as mpl
mpl.use('Agg')

import numpy as np
import xrt.runner as xrtr
import xrt.backends.raycing.screens as rsc
import xrt.plotter as xrtp
import xrt.backends.raycing.run as rr
import xrt.backends.raycing.materials as rm
import xrt.backends.raycing as raycing

from lenses import polycapillary as pl
from sources import special_sources as ss

# Create source
# Shine on a capillary
# show on screen
# done!

# Constant materials FIXME - wrap this up
mGold   = rm.Material('Au', rho=19.3)
mGlass  = rm.Material(('Si', 'O'), quantities=(1, 2), rho=2.2)

class StraightCapillaryTest(object):
    """ Implement shining through a single pipe-like object """
    def __init__(self):
        """ Tania przestrzen reklamowa """

        # Default capillary position and radius
        # must be held by the test-object
        self.x_entrance = 0.0
        self.z_entrance = 0.0
        self.R_in = 0.5

        # Default capillary dimensions: (constant in this test)
        # 10cm capillary starting at 40mm
        self.y_entrance = 40
        self.y_outrance = 140

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
                                     (0, self.y_outrance, 0))

        self.beamLine.farScreen = rsc.Screen(self.beamLine,
                                     'FarScreen',
                                     (0, self.far_screen_dist, 0))

    def make_capillary(self):
        """ Single capillary construction """
        self.capillary = pl.StraightCapillary(self.beamLine,
                                      'straightcap',
                                      x_entrance = self.x_entrance,
                                      z_entrance = self.z_entrance,
                                      y_entrance = self.y_entrance,
                                      y_outrance = self.y_outrance,
                                      R_in = self.R_in)

    def make_source(self):
        """ Many source parameters are set here by hand """

        # Source parameters
        nrays       = 2000
        E0          = 9000
        # x-direction
        distx       = 'flat'
        dx          = self.R_in * 1.5
        distxprime  = 'flat'
        dxprime     = 0.02
        # z-direction
        distz       = 'flat'
        dz          = self.R_in * 1.5
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
        """ Away from outrance """
        self.far_screen_dist = self.y_outrance + dist

    def set_capillary_entrance(self, x, z):
        """ """
        self.x_entrance = x
        self.z_entrance = z

    def set_capillary_radius(self, rad):
        """ This is necessary """
        self.R_in = rad

    def set_capillary_length(self, val):
        """ Can't change the startpoint """
        if val is 0:
            print 'Capillary length too short'
            return
        self.y_outrance = self.y_entrance + val

    def set_prefix(self, fix):
        """ Przestrzen reklamowa """
        self.prefix = str(fix)
        print 'new prefix: ', fix

    def make_run_process(self):
        """ Overloads xrt method for photon generation """
        def local_process(beamLine, shineOnly1stSource=False):
            beamTotal = None
            x = self.capillary.entrance_x()
            z = self.capillary.entrance_z()
            hitpoint = (x, 40, z)
            beamSource = beamLine.sources[0].shine(hitpoint=hitpoint)

            beamTotal, _ = self.capillary.multiple_reflect(beamSource,\
                                              maxReflections=50)

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

        # TODO - wrap this into:
        # plots = self.prepare_plots()

        # Screen at the end of capillary
        limits_near_x = [self.x_entrance - 1.2 * self.R_in,
                         self.x_entrance + 1.2 * self.R_in]
        limits_near_z = [self.z_entrance - 1.2 * self.R_in,
                         self.z_entrance + 1.2 * self.R_in]
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
                               limits=[0,20])
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
                               limits=[0,20])
            ) # different plot ends here
        plot.title = 'Divergence observations'
        plot.saveName = 'png/tests/far/' + self.prefix + '_far.png'
        plots.append(plot)

        xrtr.run_ray_tracing(plots, repeats=20, beamLine=self.beamLine,\
                processes=1)

class TaperedCapillaryTest(StraightCapillaryTest):
    """ This class is supposed to help with testing more complex
    capillary shapes: with straight axis and varying radius """

    def __init__(self):
        """ Tania przestrzen reklamowa """

        StraightCapillaryTest.__init__(self)

        # Initialize to a straight capillary
        self.R_out = self.R_in

    def set_rin_rout(self, rin, rout):
        """ riro setter """
        self.R_in = rin
        self.R_out = rout

    def make_capillary(self):
        """ Single capillary construction """
        self.capillary = pl.LinearlyTapered(self.beamLine,
                                      'straightcap',
                                      x_entrance = self.x_entrance,
                                      z_entrance = self.z_entrance,
                                      y_entrance = self.y_entrance,
                                      y_outrance = self.y_outrance,
                                      R_in = self.R_in,
                                      R_out = self.R_out)

# Similar function for each type of test would be great

def test_straight():
    """ Full test """
    # test = StraightCapillaryTest()
    test = TaperedCapillaryTest()

    # Define values to be tested
    radiuses  = [0.1, 0.2, 0.5, 0.7, 1.0]
    positions = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    lengths   = [10, 30, 50, 80, 100, 120, 150, 200]

    radiuses = [0.5]
    x_positions = [0.0, 0.3, 0.5, -0.4, 0.2, -0.3, 0.8]
    z_positions = [0.0, -0.1, 0.3, 0.2, -0.4]
    lengths = [100]
    screen_dsits = [40]

    # Run
    for radius in radiuses:
        test.set_capillary_radius(radius)
        for x_position in x_positions:
            for z_position in z_positions:
                test.set_capillary_entrance(x_position, z_position)
                for length in lengths:
                    test.set_capillary_length(length)
                    for dist in screen_dsits:
                        test.set_far_screen_distance(dist)
                        fix = 'rad_' + str(radius)
                        fix += '___pos_x_{0}_z_{1}_'.format(x_position,
                                                            z_position)
                        fix += '___len_' + str(1000+length)
                        fix += '___fsd_' + str(1000+dist)
                        test.set_prefix(fix)
                        test.run_it()

def test_tapered():
    """ Full test of a tapered capillary class """
    test = TaperedCapillaryTest()

    # Define values to be tested
    radiuses = [0.5]
    x_positions = [0.0]
    z_positions = [0.0]
    lengths = [100]
    screen_dsits = [40]
    radius_ratios = [0.5, 0.8, 0.9, 1.0, 1.1, 1.2, 1.5]

    # Run
    for radius in radiuses:
        test.set_capillary_radius(radius)
        for x_position in x_positions:
            for z_position in z_positions:
                test.set_capillary_entrance(x_position, z_position)
                for length in lengths:
                    test.set_capillary_length(length)
                    for dist in screen_dsits:
                        test.set_far_screen_distance(dist)
                        for ratio in radius_ratios:
                            test.set_rin_rout(radius, ratio * radius)
                            fix = 'Rin_' + str(radius)
                            fix += '___Rout_{0}'.format(ratio * radius)
                            fix += '___pos_x_{0}_z_{1}_'.format(x_position,
                                                                z_position)
                            fix += '___len_' + str(1000+length)
                            fix += '___fsd_' + str(1000+dist)
                            test.set_prefix(fix)
                            test.run_it()