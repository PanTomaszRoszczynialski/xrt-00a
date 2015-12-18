import xrt.backends.raycing.oes as roe
import numpy as np


class Capillary(roe.OE):
    """ Single light transmitting pipe """
    def __init__(self, *args, **kwargs):
        """ Abstract init """

        # Get y-dimension limits - by default capillary 
        # stretches from 10 to 100 mm
        self.y_entrance = kwargs.pop('y_entrance', 10)
        self.y_outrance = kwargs.pop('y_outrance', 100)

        # Update kwargs to fit xrt.OpticalElement constructor
        kwargs.update({'limPhysY' : (self.y_entrance, self.y_outrance)})

        # Get other entrance cooridinatesn
        self.x_entrance = kwargs.pop('x_entrance', 0)
        self.z_entrance = kwargs.pop('z_entrance', 0)
        # x and z outrance are defined by the capillary shape
        # and should not be of any concern here

        # We also need entrance coordinates
        # in the polar perspective
        self.phi_entrance = np.arctan2(self.x_entrance,
                                       self.z_entrance)-np.pi/2
        self.r_entrance   = np.sqrt(self.x_entrance**2 + self.z_entrance**2)

        # In xrt phi_entrance is equal to the roll
        kwargs.update({'roll' : self.phi_entrance})

        # Capital letter R will be denoting capillary 
        # radius rather then polar coordinate
        self.R_in = kwargs.pop('R_in', 1)

        # Init parent class 
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
        # FIXME - i still am not sure what is all of this doing
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

    def entrance_x(self):
        """ Basic getter """
        return self.x_entrance

    def entrance_z(self):
        """ Basic getter """
        return self.z_entrance

    def entrance_point(self):
        """ Returns cartesian coordinates of element's position """
        return self.x_in * np.sin(self.phi), self.x_in * np.cos(self.phi)

    def start_pos(self):
        """ Returns y-distance from the origin to the beginning """
        return self.y_entrance

    def end_pos(self):
        """ Returns y-distance from the origin to the finish """
        return self.y_outrance

class StraightCapillary(Capillary):
    """ Implements straight capillary parallel to the beam """
    def __init__(self, *args, **kwargs):

        # Init parent capillary class
        Capillary.__init__(self, *args, **kwargs)

        # Set straight line factors of 5 powers 
        # Separation of bent/straight capillaries
        # could be cleaner TODO
        self.p = [self.r_entrance, 0, 0, 0, 0, 0]

        # Same concept with radius
        self.pr = [self.R_in, 0, 0]

class LinearlyTapered(Capillary):
    """ You have to set the radius coefficients yourself """
    def __init__(self, *args, **kwargs):
        self.phi = kwargs['roll']
        self.x_in = kwargs.pop('x_in')
        self.p = [self.x_in, 0, 0, 0, 0, 0]

        # As a default set capillary with outward radius
        # half of the inward radius
        self.set_rin_rout(self.rin, self.rin/2.0)

        # Init parent capillary class
        Capillary.__init__(self, *args, **kwargs)

    def set_rin_rout(self, rin, rout):
        """ Set gradient to the radius """
        self.r_in = rin
        self.r_out = rout

        # Linear equation giving radius of position
        a = 1.0*(self.r_out - self.r_in) / (self.end_pos() - self.start_pos())
        b = self.r_out - self.end_pos() * a;
        self.pr = [b, a, 0]
