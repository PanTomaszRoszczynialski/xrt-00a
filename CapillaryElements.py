import numpy as np
import xrt.backends.raycing.oes as roe


# Create pythonic dictionaries to store struct like data 
# about separate sections of logic, like
# radius = {'in' : 0.1,
class BentCapillary(roe.OE):
    def __init__(self, *args, **kwargs):
        self.p          = kwargs.pop('curveCoeffs')
        self.radius     = kwargs.pop('radius')
        # Init parent Optical Element class
        roe.OE.__init__(self, *args, **kwargs)
#        self.
        kwargs.


class PolyCapillaryLens(object):
    def __init__(self, **kwargs):
        self.
