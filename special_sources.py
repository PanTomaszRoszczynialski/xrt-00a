import xrt.backends.raycing as raycing
import xrt.backends.raycing.sources as rs
from xrt.backends.raycing.sources import *
import numpy as np

defaultEnergy = 9.0e3

class DirectedSource(object):
    """Implements a geometric source with given direction"""

    def __init__(
        self, bl, name, center=(0, 0, 0), nrays=raycing.nrays, distx='normal',
        dx=0.32, disty=None, dy=0, distz='normal', dz=0.018,
        distxprime='normal', dxprime=1e-3, distzprime='normal', dzprime=1e-4,
        distE='lines', energies=(defaultEnergy,),
            polarization='horizontal', filamentBeam=False):

        self.bl = bl
        bl.sources.append(self)
        self.ordinalNum = len(bl.sources)
        self.name = name
        self.center = center  # 3D point in global system
        self.nrays = nrays

        self.distx = distx
        self.dx = dx
        self.disty = disty
        self.dy = dy
        self.distz = distz
        self.dz = dz
        self.distxprime = distxprime
        self.dxprime = dxprime
        self.distzprime = distzprime
        self.dzprime = dzprime
        self.distE = distE
        if self.distE == 'lines':
            self.energies = np.array(energies)
        else:
            self.energies = energies
        self.polarization = polarization
        self.filamentBeam = filamentBeam

    def _apply_distribution(self, axis, distaxis, daxis):
        if (distaxis == 'normal') and (daxis > 0):
            axis[:] = np.random.normal(0, daxis, self.nrays)
        elif (distaxis == 'flat'):
            if raycing.is_sequence(daxis):
                aMin, aMax = daxis[0], daxis[1]
            else:
                if daxis <= 0:
                    return
                aMin, aMax = -daxis*0.5, daxis*0.5
            axis[:] = np.random.uniform(aMin, aMax, self.nrays)
#        else:
#            axis[:] = 0

    def shine(self, hitpoint=(0,10,0), toGlobal=True,\
              withAmplitudes=False, accuBeam=None):
        """Returns the source beam. If *toGlobal* is True, the output is in
        the global system. If *withAmplitudes* is True, the resulted beam
        contains arrays Es and Ep with the *s* and *p* components of the
        electric field.
        """
        bo = Beam(self.nrays, withAmplitudes=withAmplitudes)  # beam-out
        bo.state[:] = 1
# =0: ignored, =1: good,
# =2: reflected outside of working area, =3: transmitted without intersection
# =-NN: lost (absorbed) at OE#NN (OE numbering starts from 1!) If NN>1000 then
# the slit with ordinal number NN-1000 is meant.

# in local coordinate system:
        self._apply_distribution(bo.y, self.disty, self.dy)

        self._apply_distribution(bo.x, self.distx, self.dx)
        self._apply_distribution(bo.z, self.distz, self.dz)

        # Only momentum distribution is different then in normal
        # source
        normfactor = np.sqrt(sum([x**2 for x in hitpoint]))
        a0 = 1. * hitpoint[0] / normfactor
        aMin, aMax = a0 - self.dxprime, a0 + self.dxprime
        bo.a = np.random.uniform(aMin, aMax, self.nrays)

        b0 = 1. * hitpoint[2] / normfactor
        bMin, bMax = b0 - self.dzprime, b0 + self.dzprime
        bo.c = np.random.uniform(bMin, bMax, self.nrays)

# normalize (a,b,c):
        ac = bo.a**2 + bo.c**2
        if sum(ac > 1) > 0:
            bo.b[:] = (ac + 1)**0.5
            bo.a[:] /= bo.b
            bo.c[:] /= bo.b
            bo.b[:] = 1.0 / bo.b
        else:
            bo.b[:] = (1 - ac)**0.5
        if self.distE is not None:
            if accuBeam is None:
                bo.E[:] = make_energy(self.distE, self.energies, self.nrays,
                                      self.filamentBeam)
            else:
                bo.E[:] = accuBeam.E[:]
        make_polarization(self.polarization, bo, self.nrays)
        if toGlobal:  # in global coordinate system:
            raycing.virgin_local_to_global(self.bl, bo, self.center)
        return bo
