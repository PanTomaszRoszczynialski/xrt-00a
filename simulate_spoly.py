"""
First try to recreate matlab code generating
lens' hexagonal structure
"""
import numpy as np
import matplotlib.pyplot as plt
from random import random as rand

class HexStructure(object):
    def __init__(self, capillary_diameter = 0.0025,\
                nx_capillary = 9,\
                ny_bundle = 7):
        # Outer diameter (touching)
        self.capillary_diameter = capillary_diameter

        # Number of bundles along vertical direction (must be odd)
        self.ny_bundle = ny_bundle

        # Number of capillaries in bundle along 
        # horizontal direction (must be odd)
        self.nx_capillary = nx_capillary

        # Check oddity
        if self.ny_bundle % 2 == 0:
            print "Number of bundles must be odd, adding one"
            self.ny_bundle += 1
        if self.nx_capillary %2 ==0:
            print "Number of capillaries must be odd, adding one"
            self.nx_capillary += 1

        # Inner diameter (x-ray microsource ? )
        self.channel_diameter = self.capillary_diameter * 0.5

        # Spacing between capilaries (touching)
        self.bundlespacing = (self.nx_capillary - 1) *\
                                self.capillary_diameter *\
                                np.sqrt(3)/2 + self.capillary_diameter

        self.capillary_lens_xy()

    def capillary_bundle_xy(self, xbundle, ybundle, sigma):
        index = 0
        nxpol_capillary = (self.nx_capillary - 1)/2

        atemp = self.capillary_diameter * nxpol_capillary

        # Prepare returned lists
        xci = []
        yci = []
        for ix in range(-2*nxpol_capillary, 2*nxpol_capillary +1):
            for iy in range(-2*nxpol_capillary, 2*nxpol_capillary +1):

                x0 = self.capillary_diameter * ix +\
                    self.capillary_diameter/2*iy + \
                    sigma * (rand() - 0.5) * \
                    (self.capillary_diameter - self.channel_diameter)

                y0 = np.sqrt(3)/2 * self.capillary_diameter * iy +\
                    sigma * (rand() - 0.5) * \
                    (self.capillary_diameter - self.channel_diameter)

                in_bundle = self.isInHexagon(y0,x0,\
                                self.capillary_diameter*nxpol_capillary)

                if in_bundle:
                    xci.append(xbundle + x0)
                    yci.append(ybundle + y0)
        return xci, yci



    def capillary_lens_xy(self, sigma_position = 0.1):

        # This has to be an integer, do not add dots
        nypol_bundle = (self.ny_bundle - 1)/2

        # Prepare return containers
        xci = []
        yci = []
        xi = []
        yi = []

        # Add + 1 because range is 0 based...
        for ix in range(-2*nypol_bundle, 2*nypol_bundle +1):
#            print str(ix)
            for iy in range(-2*nypol_bundle, 2*nypol_bundle +1):

                x0 = np.sqrt(3)/2.0 * self.bundlespacing * iy +\
                        0 * sigma_position * (rand()-0.5) *\
                        (self.capillary_diameter - self.channel_diameter)

                y0 = self.bundlespacing * ix + \
                        self.bundlespacing * iy / 2.0 + \
                        0 * sigma_position * (rand()-0.5) *\
                        (self.capillary_diameter - self.channel_diameter)

                in_lens = self.isInHexagon(y0, x0, self.bundlespacing *\
                                            nypol_bundle)

                if in_lens:
                    #print x0, y0
                    xci0, yci0 = self.capillary_bundle_xy(x0, y0,\
                                    sigma_position)
                    xci.append(xci0)
                    yci.append(yci0)
                    xi.append(x0)
                    yi.append(y0)

        self.xi = xi
        self.yi = yi
        self.xci = xci
        self.yci = yci
        return xi, yi, xci, yci


    def isInHexagon(self, y0, x0, d):
        tol = 1.001

        war1 = abs(y0) <= d * tol * np.sqrt(3)/2
        war2 = abs(np.sqrt(3)/2* x0 + 1/2. * y0) <= tol * d * np.sqrt(3)/2
        war3 = abs(np.sqrt(3)/2* x0 - 1/2. * y0) <= tol * d * np.sqrt(3)/2

        return war1 and war2 and war3

    def test(self):
        print len(self.xci)
        plt.plot(self.xi, self.yi,'ko')
        plt.show()


if __name__ == '__main__':

    hello = HexStructure(nx_capillary=14, ny_bundle=13 )
    hello.test()
