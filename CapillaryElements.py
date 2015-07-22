import numpy as np
import xrt.backends.raycing.oes as roe


# Create pythonic dictionaries to store struct like data 
# about separate sections of logic, like
# radius = {'in' : 0.1}
class BentCapillary(roe.OE):
    def __init__(self, *args, **kwargs):
        # Polynomial factors defining shape in the y direction
        self.p          = kwargs.pop('curveCoeffs')

        # Factors defining local radius of the capillary
        # rIn, rOut, rMax
        self.radius     = kwargs.pop('radius')

        # Parameters set for the whole Lens

        # Init parent class (Optical Element )
        roe.OE.__init__(self, *args, **kwargs)


class PolyCapillaryLens(object):
    def __init__(self, **kwargs):
        # Hopefuly it's save to make a beamLine member
        self.beamLine   = kwargs.pop('beamLine')
        # Elements' positions in y direction
        self.y          = kwargs.pop('y_setting')
        # Lens defining diameters (in, out, max)
        self.D          = kwargs.pop('D_setting')
        # Material of capillaries
        self.material   = kwargs.pop('material', None)

    def setStructure(self, structure):
        # Structure defining parameters in 'x' and 'z' coordinates
        # No loading options needed: it will be done outside this
        # class at instantiation 
        self.structure = structure

    def capillaryParameters(self, x_in, roll):
        args = [self.beamLine, 'bent', [0,0,0]]
        # Position of capillary in polar coordinates
        kwargs = {'roll' : roll, 'x_in' : x_in}

        # Physical limit in y direction
        limPhysY = [self.y['y1'], self.y['y2'] ]
        kwargs.update({'limPhysY' : limPhysY})

        # Radius settings
        rIn = self.structure.getCapillaryRadius()
        rOut = rIn * self.D['Dout'] / self.D['Din']
        rMax = rIn * self.D['Dmax'] / self.D['Din']
        kwargs.update({'rIn' : rIn, 'rOut' : rOut, 'rMax' : rMax})

        # Parameters needed for capillary shape in z direction
        kwargs.update(self.y)
        kwargs.update(self.D)

        # Material of capillary
        kwargs.update({'material' : self.material})

        return args, kwargs

    def setCapillarySettings(self, settings):
        self.c_radius = settings['radius']
        self.c_wall = settings['wall']

    def getCapillaries(self, beamLine):
        # Prepare containers
        capillaries = []
        toPlot = []

        # Generate capillaries with positions given
        # by the lens structure
        for r, phi in structure.genPolars():
            roll = phi
            x    = r

            # Capillary should care only about x and phi variable
            args, kwargs = capillaryParameters(x, roll)
            capillary = BentCapillary(*args, **kwargs)
            capillaries.append(capillary)

        beamline.capillaries = capillaries
        beamLine.toPlot = toPlot

        return beamLine

