import xrt.backends.raycing.oes as roe
import xrt.backends.raycing.materials as rm
import numpy as np

# Constant materials
mGold   = rm.Material('Au', rho=19.3)
mGlass  = rm.Material(('Si', 'O'), quantities=(1, 2), rho=2.2)

class Capillary(roe.OE):
    """ Single light transmitting pipe """
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
        # FIXME: this method seems to be generating highly
        # erroneus behavior -- photons get outside of the
        # capillary and also bounce more times then allowed
        # in the xrt.multiple_reflect() method
        den = np.cos(np.arctan(self.local_x0Prime(s)))**2
        return self.local_r0(s) / (np.cos(phi)**2/den + np.sin(phi)**2)

    def local_n(self, s, phi):
        a = -np.sin(phi)
        b = -np.sin(phi)*self.local_x0Prime(s) - self.local_r0Prime(s)
        c = -np.cos(phi)
        norm = np.sqrt(a**2 + b**2 + c**2)
        return a/norm, b/norm, c/norm

    def xyz_to_param(self, x, y, z):
        """ *s*, *r*, *phi* are cylindrc-like coordinates of the capillary.
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

    def entrance_point(self):
        """ Returns cartesian coordinates of element's position """
        return self.x_in * np.sin(self.phi), self.x_in * np.cos(self.phi)

class StraightCapillary(Capillary):
    """ Implements straight capillary parallel to the beam """
    def __init__(self, *args, **kwargs):
        self.x_in   = kwargs.pop('x_in')
        self.r_in   = kwargs.pop('r')

        # This should be present only when y-information
        # are used for shape coefficients calculations
        y = kwargs['limPhysY']
        self.y0 = y[0]
        self.y1 = y[1]
        self.phi = kwargs['roll']

        # Set straight line factors of 5 powers 
        # Separation of bent/straight capillaries
        # could be cleaner TODO
        self.p = [self.x_in, 0, 0, 0, 0, 0]

        # Same concept with radius
        self.pr = [self.r_in, 0, 0]

        # This should be settable from the testing module FIXME
        # kwargs.update({'material' : mGold})
        kwargs.update({'material' : mGlass})

        # Init parent capillary class
        Capillary.__init__(self, *args, **kwargs)

class LinearlyTapered(Capillary):
    """ You have to set the radius coefficients yourself """
    def __init__(self, *args, **kwargs):
        self.phi = kwargs['roll']
        self.x_in = kwargs.pop('x_in')
        self.p = [self.x_in, 0, 0, 0, 0, 0]

        # As a deafault this is a Straight capillary
        self.pr = [self.r_in, 0, 0]

        # This should be settable from the testing module FIXME
        # kwargs.update({'material' : mGold})
        kwargs.update({'material' : mGlass})

        # Init parent capillary class
        Capillary.__init__(self, *args, **kwargs)

    def set_rin_rout(self, rin, rout):
        """ Set gradient to the radius """
        self.rin = rin
        self.rout = rout
