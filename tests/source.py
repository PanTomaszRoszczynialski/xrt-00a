import numpy as np
from sources import special_sources as ss
import xrt.runner as xrtr
import xrt.backends.raycing.screens as rsc
import xrt.plotter as xrtp
import xrt.backends.raycing.run as rr
import xrt.backends.raycing as raycing

# TODO Move those 
rIn = 0.5

class SourceTest(object):
    """ Class for testing sources (duh..) """
    def __init__(self, source):
        """ Takes pointer to source constructor """
        # xrt commons
        self.beamLine = raycing.BeamLine(height=0)

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

            beamTotal = None

            # This structure can be re-set
            for hp in self.hitpoints:

                # And this should read hp = structure.next()
                hitpoint = (hp[0], 40, hp[1])

                # Generate photons
                beamSource = beamLine.sources[0].shine(hitpoint=hitpoint)

                # Accumulate
                if beamTotal is None:
                    beamTotal = beamSource
                else:
                    beamTotal.concatenate(beamSource)

            # Expose screens
            sourceScreen = beamLine.entScreen.expose(beamTotal)
            farScreen = beamLine.farScreen.expose(beamTotal)

            # Show beamlines after exposition
            outDict = {'SourceScreen' : sourceScreen}
            outDict['FarScreen'] = farScreen

            return outDict

        rr.run_process = source_process

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

        xrtr.run_ray_tracing(plots, repeats=100, beamLine=beamLine,\
                processes=1)

def run_test():
    """ Creates circular arangement of illuminated regions
    and plots xrt histograms right at the source end,
    and also somewhere far (to observe beam divergence) """
    # source = ss.FittedSource
    source = ss.DirectedSource
    test = SourceTest(source)

    # Set circular structer, as a test
    structure = []
    for it in range(-20, 20):
        r = 15
        phi = 2. * np.pi * it / 40.0
        x = r * np.cos(phi)
        y = r * np.sin(phi)
        structure.append([x, y])
    test.set_structure(structure)
    test.make_run_process()
    test.run_test()
