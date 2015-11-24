import numpy as np
from sources import special_sources as ss
import xrt.runner as xrtr
import xrt.backends.raycing.screens as rsc
import xrt.plotter as xrtp
import xrt.backends.raycing.run as rr
import xrt.backends.raycing as raycing

# TODO reorganize more
bins = 256
rIn = 0.5

class SourceTest(object):
    def __init__(self, source):
        """ Take pointer to source constructor """
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

        # Create source
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

    def make_run_process(self):
        """ Overload xrt method managing screens """
        def fitted_source_process(beamLine, shineOnly1stSource=False):

            beamTotal = None

            for it in range(-20,20):
                r = 15
                phi = 2.*np.pi * it / 40.0
                x = r * np.cos(phi)
                y = r * np.sin(phi)
                hitpoint = (x, 40, y)
                beamSource = beamLine.sources[0].shine(hitpoint=hitpoint)

                if beamTotal is None:
                    beamTotal = beamSource
                else:
                    beamTotal.concatenate(beamSource)

            sourceScreen = beamLine.entScreen.expose(beamTotal)
            farScreen = beamLine.farScreen.expose(beamTotal)

            outDict = {'SourceScreen' : sourceScreen}
            outDict['FarScreen'] = farScreen

            return outDict

        rr.run_process = fitted_source_process

    def run_test(self):
        beamLine = self.beamLine

        limits = [-20,20]
        plots = []

        # Screen at the source beg
        plot = xrtp.XYCPlot(
            # Using named parameters might be good here TODO
            'SourceScreen', (1, 3,),
            # This is still inside plot
            xaxis=xrtp.XYCAxis(r'$x$',
            'mm', bins=bins, ppb=2, limits=limits),
            # As is this
            yaxis=xrtp.XYCAxis(r'$z$', 'mm', bins=bins,
            ppb=2, limits=limits),
            caxis='category', beamState='SourceScreen')
        plot.title = 'Source structure'
        plots.append(plot)

        plot = xrtp.XYCPlot(
            'FarScreen', (1, 3,),
            xaxis=xrtp.XYCAxis(r'$x$', 'mm', bins=bins, ppb=2, limits=None),
            yaxis=xrtp.XYCAxis(r'$z$', 'mm', bins=bins, ppb=2, limits=None),
            caxis='category', beamState='FarScreen')
        plot.title = 'Source divergance'
        plots.append(plot)
        xrtr.run_ray_tracing(plots, repeats=100, beamLine=beamLine,\
                processes=1)

def test_fitted_source():
    beamLine = raycing.BeamLine(height=0)
    ss.FittedSource(
        beamLine,'DirectedSource',(0,40,0), nrays=nrays,
        distx=distx, dx=dx, distxprime=distxprime, dxprime=dxprime,
        distz=distz, dz=dz, distzprime=distzprime, dzprime=dzprime,
        distE='lines', energies=(E0,), polarization='horizontal')

    # Insert screen at the screen position
    beamLine.entScreen = rsc.Screen(beamLine, 'EntranceScreen',(0,40,0))
    beamLine.farScreen = rsc.Screen(beamLine, 'FarScreen', (0,11120,0))

    return beamLine

def fitted_source_process(beamLine, shineOnly1stSource=False):
    beamTotal = None

    for it in range(-20,20):
        r = 15
        phi = 2.*np.pi * it / 40.0
        x = r * np.cos(phi)
        y = r * np.sin(phi)
        hitpoint = (x, 40, y)
        beamSource = beamLine.sources[0].shine(hitpoint=hitpoint)

        if beamTotal is None:
            beamTotal = beamSource
        else:
            beamTotal.concatenate(beamSource)

    sourceScreen = beamLine.entScreen.expose(beamTotal)
    farScreen = beamLine.farScreen.expose(beamTotal)

    outDict = {'SourceScreen' : sourceScreen}
    outDict['FarScreen'] = farScreen

    return outDict


def run_test():
    source = ss.FittedSource
    test = SourceTest(source)
    test.make_run_process()
    test.run_test()
