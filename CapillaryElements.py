import numpy as np
from LensPolynomial import getPolyCoeffs, getRadiusCoeffs
import xrt.backends.raycing.oes as roe
import xrt.backends.raycing.materials as rm

# Constant materials
mGold   = rm.Material('Au', rho=19.3)

# Creating pythonic dictionaries to store struct like data 
# about separate sections of logic, like
# radius = {'in' : 0.1}
class Capillary(roe.OE):
    def __init__(self, *args, **kwargs):
        # Init parent class (Optical Element )
        roe.OE.__init__(self, *args, **kwargs)
        self.isParametric = True

    def local_x0(self, s):
        return self.p[0] + self.p[1]*s + self.p[2]*s**2\
                + self.p[3]*s**3 + self.p[4]*s**4 + self.p[5]*s**5

    def local_x0Prime(self, s):
        return self.p[1] + 2*self.p[2]*s + 3*self.p[3]*s**2\
                + 4*self.p[4]*s**3 + 5*self.p[5]*s**4

    def local_r0(self, s):
        return self.pr[0] + self.pr[1]*s + self.pr[2]*s**2

    def local_r0Prime(self,s):
        return self.pr[1] + 2*self.pr[2]*s

    def local_r(self, s, phi):
        den = np.cos(np.arctan(self.local_x0Prime(s)))**2
        return self.local_r0(s) / (np.cos(phi)**2/den + np.sin(phi)**2)

    def local_n(self, s, phi):
        a = -np.sin(phi)
        b = -np.sin(phi)*self.local_x0Prime(s) - self.local_r0Prime(s)
        c = -np.cos(phi)
        norm = np.sqrt(a**2 + b**2 + c**2)
        # FIXME: a and c probably should also get minues
        # but due to symmetry it's not visible for now?
        return a/norm, -b/norm, c/norm

    def xyz_to_param(self, x, y, z):
        """ *s*, *r*, *phi* are cynindrc-like coordinates of the capillary.
        *s* is along y in inverse direction, started at the exit,
        *r* is measured from the capillary axis x0(s)
        *phi* is the polar angle measured from the z (vertical) direction."""
        s = y
        phi = np.arctan2(x - self.local_x0(s), z)
        r = np.sqrt((x-self.local_x0(s))**2 + z**2)
        return s, phi, r

    def param_to_xyz(self, s, phi, r):
        x = self.local_x0(s) + r*np.sin(phi)
        y = s
        z = r * np.cos(phi)
        return x, y, z

class BentCapillary(Capillary):
    def __init__(self, *args, **kwargs):
        # Prepare z direction polynomial of focusing capillary
        y = kwargs.pop('y')
        D = kwargs.pop('D')
        x_in = kwargs.pop('x_in')
        self.p = getPolyCoeffs(x_in, y, D)

        # Prepare variable radius
        r = kwargs.pop('radius')
        self.pr = getRadiusCoeffs(y, r)

        # Init parent capillary class
        Capillary.__init__(self, *args, **kwargs)


class StraightCapillary(Capillary):
    def __init__(self, *args, **kwargs):
        self.x_in   = kwargs.pop('x_in')
        self.r_in   = kwargs.pop('r')
        # No need to repeat values at instantiation
        y = kwargs['limPhysY']
        self.y0 = y[0]
        self.y1 = y[1]

        # Set straight line factors of 5 powers 
        # Separation of bent/straight capillaries
        # could be cleaner TODO
        self.p = [self.x_in, 0, 0, 0, 0, 0]

        # Same concept with radius
        self.pr = [self.r_in, 0, 0]

        # Init parent capillary class
        Capillary.__init__(self, *args, **kwargs)

class Pinhole(StraightCapillary):
    def __init__(self, *args, **kwargs):
        self.pinlen = 0.01 #const
        self.y_in   = kwargs.pop('y_in')
        limPhysY = [self.y_in - self.pinlen, self.y_in]
        kwargs.update({'limPhysY' : limPhysY})
        kwargs.update({'material' : mGold})

        # Init parent classes
        StraightCapillary.__init__(self, *args, **kwargs)

class PolyCapillaryLens(object):
    def __init__(self, **kwargs):
        # Elements' positions in y direction
        self.y          = kwargs.pop('y_settings')
        # Lens defining diameters (in, out, max)
        self.D          = kwargs.pop('D_settings')
        # Material of capillaries
        self.material   = kwargs.pop('material', None)

    def setStructure(self, structure):
        # Structure defining parameters in 'x' and 'z' coordinates
        # No loading options needed: it will be done outside this
        # class at instantiation 
        self.structure = structure

    def capillaryParameters(self, x_in, roll, beamLine):
        # Default 'positional' parameters
        args = [beamLine, 'bent', [0,0,0]]

        # Named parameters (dict of them)
        # Position of capillary in polar coordinates
        kwargs = {'roll' : roll, 'x_in' : x_in}

        # Physical limit in y direction
        limPhysY = [self.y['y1'], self.y['y2'] ]
        kwargs.update({'limPhysY' : limPhysY})

        # Radius settings
        rIn = self.structure.getCapillaryRadius()
        rOut = rIn * self.D['Dout'] / self.D['Din']
        rMax = rIn * self.D['Dmax'] / self.D['Din']
        radius = {'rIn' : rIn, 'rOut' : rOut, 'rMax' : rMax}
        kwargs.update({'radius' : radius})

        # Parameters needed for capillary shape in z direction
        kwargs.update({'y' : self.y})
        kwargs.update({'D' : self.D})

        # Material of capillary
        kwargs.update({'material' : self.material})

        return args, kwargs

    def getCapillaries(self, beamLine):
        # Prepare containers
        capillaries = []
        toPlot = []

        # Generate capillaries with positions given
        # by the lens structure
        for r, phi in self.structure.genPolars():
            roll = phi
            x    = r

            # Capillary should care only about x and phi variable
            args, kwargs = self.capillaryParameters(x, roll, beamLine)
            capillary = BentCapillary(*args, **kwargs)
            capillaries.append(capillary)

            # Save capillaries shown on z=0 coss-section ? 
            # Z = 0 is no longer special
            # and as it is clear neither is phi = pi/3,
            # so some other idea for crosssection plot
            # is needed TODO
            # DEBUG quick polar to cartesian re-transformation
            x_cap = r * np.cos(phi)
            if abs(x_cap) < 0.005:
                toPlot.append(len(capillaries))

        beamLine.capillaries = capillaries
        beamLine.toPlot = toPlot

        # beamLine should've been a reference to some outside 
        # object, but who knows
        return beamLine

