""" Number of processes depends on the machine,
    mpl.use('Agg') should be uncommented when You want
    to save results without displaying any plots """
# SETME
_processes = 1
import matplotlib as mpl
mpl.use('Agg')

import numpy as np
from sources import special_sources as ss
import xrt.backends.raycing.sources as rs
import xrt.runner as xrtr
import xrt.backends.raycing.screens as rsc
import xrt.plotter as xrtp
import xrt.backends.raycing.run as rr
import xrt.backends.raycing as raycing

# Almost obsolete
class SourceTest(object):
    """ Class for testing sources """
    def __init__(self, source):
        """ Takes pointer to source constructor """
        # xrt commons
        self.beamLine = raycing.BeamLine(height=0)

        # New idea
        self.beamTotal = None

        rIn = 0.5
        # Source parameters
        nrays       = 300
        E0          = 9000
        # x-direction
        distx       = 'flat'
        dx          = rIn * 1.5
        distxprime  = 'flat'
        dxprime     = 0.02
        # z-direction
        distz       = 'flat'
        dz          = rIn * 1.5
        distzprime  = 'flat'
        dzprime     = 0.02

        # Initialize tested structure with a line of source hitpoints,
        # they are used only in the run_process mechanism
        structure = []
        for it in range(-10, 10):
            r = 15
            x = r * it / 5.
            y = 0
            structure.append([x, y])
        self.hitpoints = structure

        # Create source of tested type
        source(
            self.beamLine,'DirectedSource',(0,40,0), nrays=nrays,
            distx=distx, dx=dx, distxprime=distxprime, dxprime=dxprime,
            distz=distz, dz=dz, distzprime=distzprime, dzprime=dzprime,
            distE='lines', energies=(E0,), polarization='horizontal')

        # Insert screen right behind the source ...
        self.beamLine.entScreen = rsc.Screen(self.beamLine,
                                        'EntranceScreen',
                                        (0,40,0))
        # ... and somewhere far
        self.beamLine.farScreen = rsc.Screen(self.beamLine,
                                        'FarScreen',
                                        (0,400,0))

    def set_structure(self, structure):
        """ Call this method before make_run_process,
        makes a list of [x, y] pairs with positions
        we wan't to illuminate with the source """
        self.hitpoints = structure

    def make_run_process(self):
        """ Overload xrt method managing screens """
        def source_process(beamLine, shineOnly1stSource=False):

            # beamTotal = None

            # This structure can be re-set
            for hp in self.hitpoints:

                # And this should read hp = structure.next()
                hitpoint = (hp[0], 40, hp[1])

                # Generate photons
                beamSource = beamLine.sources[0].shine(hitpoint=hitpoint)

                # Accumulate
                if self.beamTotal is None:
                    self.beamTotal = beamSource
                else:
                    self.beamTotal.concatenate(beamSource)

            # Expose screens
            sourceScreen = beamLine.entScreen.expose(self.beamTotal)
            farScreen = beamLine.farScreen.expose(self.beamTotal)

            # Show beamlines after exposition
            outDict = {'SourceScreen' : sourceScreen}
            outDict['FarScreen'] = farScreen

            return outDict

        rr.run_process = source_process

    def get_beam(self):
        """ Get beamLine object holding all photon data """
        return self.beamTotal

    def run_test(self):
        beamLine = self.beamLine

        # xrt.Plot constants
        bins = 256
        limits = [-20,20]
        plots = []

        # Screen at the source beg
        plot = xrtp.XYCPlot(
            # Using named parameters might be good here TODO
            'SourceScreen', (1, 3,),
            caxis='category', beamState='SourceScreen',
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
        plot.title = 'Source structure'
        plots.append(plot)

        plot = xrtp.XYCPlot(
            'FarScreen', (1, 3,),
            xaxis=xrtp.XYCAxis(r'$x$', 'mm',
            bins=bins, ppb=2, limits=None),
            yaxis=xrtp.XYCAxis(r'$z$', 'mm',
            bins=bins, ppb=2, limits=None),
            caxis='category', beamState='FarScreen')
        plot.title = 'Source divergance'
        plots.append(plot)

        _repeats = 20
        xrtr.run_ray_tracing(plots, repeats=_repeats, beamLine=beamLine,\
                processes=1)

# This is obsolete
def run_test():
    """ Creates circular arangement of illuminated regions
    and plots xrt histograms right at the source end,
    and also somewhere far (to observe beam divergence) """
    # source = ss.FittedSource
    # source = ss.DirectedSource
    # test = SourceTest(source)
    test = GeometricSourceTest()

    # Set circular structer, as a test
    structure = []
    for it in range(-20, 20):
        r = 15
        phi = 2. * np.pi * it / 40.0
        x = r * np.cos(phi)
        y = r * np.sin(phi)
        structure.append([x, y])
    # test.set_structure(structure)
    test.make_run_process()
    test.run_it()
    beamGlobalTotal = test.get_beam()
    return beamGlobalTotal


class GeometricSourceTest(object):
    """ Implements methods for easy observations of chosen source """
    def __init__(self):
        """ Co robi traktor u fryzjera? Warkocze """

        # Testing new functionality
        self.beamTotal = None

        # Physical end of the source
        self.y_outrance = 1

        # Distant screen
        self.far_screen_dist = 5

        # Source dimensions
        self.x_size = 1
        self.z_size = 1

        # Source divergence
        self.x_divergence = 0.1
        self.z_divergence = 0.1

        # Source energy parameters
        self.distE = 'normal'
        self.energies = (9000, 100)

        # Number of photons in one iteration of one thread
        self.nrays = 1000

        # Total Number of Photons
        self.TNoP = 5000

        # Adjust number of iterations to get desired output size
        self.repeats = self.TNoP / self.nrays

        # Other
        self.prefix = 'setme'

    def get_beam(self):
        """ Get beamLine object holding all photon data """
        return self.beamTotal

    def run_it(self):
        """ Runs the whole operation """

        self.make_it()
        xrtr.run_ray_tracing(self.plots,\
                            repeats=self.repeats,\
                            beamLine=self.beamLine,\
                            processes=_processes)

    def make_it(self):
        """ Re-constructs the object """
        self.make_source()
        self.make_screens()
        self.make_run_process()
        self.make_plots()

    def make_run_process(self):
        """ Overloads xrt method for photon generation """
        def local_process(beamLine, shineOnly1stSource=False):
            beamSource = beamLine.sources[0].shine()

            # Hold photons for export
            if self.beamTotal is None:
                self.beamTotal = beamSource
            else:
                self.beamTotal.concatenate(beamSource)

            # Expose screens
            exitScreen = beamLine.exitScreen.expose(beamSource)
            farScreen = beamLine.farScreen.expose(beamSource)

            # Show beamlines after exposition
            outDict = {'ExitScreen' : exitScreen}
            outDict['FarScreen'] = farScreen

            return outDict

        rr.run_process = local_process

    def make_plots(self):
        """ tania przestrzen reklamowa """

        # Plots' resolution
        bins = 256
        self.plots = []

        plot = xrtp.XYCPlot(
            # Using named parameters might be good here TODO
            'ExitScreen', (1, 3,),
            beamState='ExitScreen',
            # This is still inside plot
            xaxis=xrtp.XYCAxis(r'$x$', 'mm',
                               bins=bins,
                               ppb=2,
                               limits=None),
            # As is this
            yaxis=xrtp.XYCAxis(r'$z$', 'mm',
                               bins=bins,
                               ppb=2,
                               limits=None),
            # Colorbar and sidebar histogram
            caxis=xrtp.XYCAxis('Energy',
                               '[eV]',
                               data=raycing.get_energy,
                               bins=bins,
                               ppb=2,
                               limits=None)
            ) # plot ends here
        plot.title = 'Near source'
        plot.saveName = 'png/tests/near/' + self.prefix + '_near.png'
        self.plots.append(plot)

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
            caxis=xrtp.XYCAxis('Energy',
                               '[eV]',
                               data=raycing.get_energy,
                               bins=bins,
                               ppb=2,
                               limits=None)
            ) # different plot ends here
        plot.title = 'Divergence observations'
        plot.saveName = 'png/tests/far/' + self.prefix + '_far.png'
        self.plots.append(plot)

    def make_screens(self):
        """ One screen at the exit and one somewhere far """
        self.beamLine.exitScreen = rsc.Screen(self.beamLine,
                                     'Exitscreen',
                                     (0, self.y_outrance, 0))

        self.beamLine.farScreen = rsc.Screen(self.beamLine,
                                     'FarScreen',
                                     (0, self.far_screen_dist, 0))

    def make_source(self):
        """ Many source parameters are set here by hand """

        # Source parameters
        nrays       = self.nrays
        distE       = self.distE
        energies    = self.energies
        # x-direction
        distx       = 'flat'
        dx          = self.x_size
        distxprime  = 'normal'
        dxprime     = self.x_divergence
        # z-direction
        distz       = 'flat'
        dz          = self.z_size
        distzprime  = 'normal'
        dzprime     = self.z_divergence

        self.polarization = None

        # Sources appends themselves to the beamline
        # so we recreate the beamline to always use bl.sources[0]
        self.beamLine = raycing.BeamLine(height=0)
        rs.GeometricSource(
            self.beamLine,'DirectedSource',(0,0,0), nrays=nrays,
            distx=distx, dx=dx, distxprime=distxprime, dxprime=dxprime,
            distz=distz, dz=dz, distzprime=distzprime, dzprime=dzprime,
            distE=distE, energies=energies,
            polarization=self.polarization)

    def set_prefix(self, fix):
        """ This helps multiple file management """
        self.prefix = fix

    def set_xz_size(self, sizex, sizez):
        """ don't panic """
        self.x_size = sizex
        self.z_size = sizez

    def set_total_number_of_photons(self, tnop):
        """ Set how many photons You want generated """
        self.TNoP = tnop

        # Adjust number of iterations to get desired output size
        self.repeats = self.TNoP / self.nrays

    # TODO - if source parameters could be controlled
    def set_xz_divergence(self, divx, divz):
        """ Sets photons distribution sigma-parameters """
        self.x_divergence = divx
        self.z_divergence = divz

    def set_energies(self, E0, sigma):
        """ In electronovolts please """
        if sigma is not 0:
            self.distE = 'normal'
            self.energies = (E0, sigma)
        else:
            self.distE = 'lines'
            self.energies = (E0,)

def test_geometric():
    """ Runs any kind of tests on geometric source """
    test = GeometricSourceTest()

    test.set_xz_size(1,1)
    test.set_total_number_of_photons(1000)
    test.set_xz_divergence(0.01,0.01)
    z_divs = [0.1 * x for x in range(1,11)]
    z_divs = [0.01]
    E0 = 9000
    e_sigmas = [10, 20, 30, 40, 50]
    for div in z_divs:
        for sigma in e_sigmas:
            test.set_energies(E0, sigma)
            test.set_xz_divergence(0.01, div)
            fix = 'zdiv_{0}'.format(div)
            fix += '___e_sigma_{0}'.format(sigma)
            test.set_prefix(fix)
            test.run_it()

def create_geometric():
    """ Create photons from a Geometric Source """
    source = GeometricSourceTest()
    source.set_xz_size(0.1, 0.1)
    source.set_xz_divergence(0.1, 0.1)
    source.set_total_number_of_photons(1e5)
    source.run_it()

    return source.get_beam()


if __name__ == '__main__':

    test_geometric()
    # GlobalTotal = create_geometric()
