import csv
import numpy as np
import matplotlib.pyplot as plt

import xrt.backends.raycing as raycing
from special_sources import DirectedSource
from special_sources import FittedSource
from xrt.backends.raycing.screens import Screen

def add_FittedSource(beamLine):
    E0 = 20000
    # Source parameters
    distx       = 'flat'
    dx          = 0.006
    distxprime  = 'flat'
    dxprime     = 0.03
    # z-direction
    distz       = 'flat'
    dz          = 0.006
    distzprime  = 'flat'
    dzprime     = 0.03
    # [0] - Source of light
    FittedSource(
        beamLine,'FittedSource',(0,40,0), nrays=3000,
        distx=distx, dx=dx, distxprime=distxprime, dxprime=dxprime,
        distz=distz, dz=dz, distzprime=distzprime, dzprime=dzprime,
        distE='lines', energies=(E0,), polarization=None)


def add_DirectedSource(beamLine):
    E0 = 20000
    # Source parameters
    distx       = 'flat'
    dx          = 1e-9
    distxprime  = 'flat'
    dxprime     = 0.002
    # z-direction
    distz       = 'flat'
    dz          = 1e-9
    distzprime  = 'flat'
    dzprime     = 0.002
    # [0] - Source of light
    DirectedSource(
        beamLine,'DirectedSource',(0,0,0), nrays=3000,
        distx=distx, dx=dx, distxprime=distxprime, dxprime=dxprime,
        distz=distz, dz=dz, distzprime=distzprime, dzprime=dzprime,
        distE='lines', energies=(E0,), polarization=None)


def view_source():
    # Make beamline
    beamLine = raycing.BeamLine(height=0)
    # [0] - Source of light
    add_FittedSource(beamLine)
    # Screen
    screen = Screen(beamLine, 'dupa', (0,40,0))
    # Multiple hitpoints
    shining = beamLine.sources[0].shine(hitpoint = (3,40,1))
    shining_b = beamLine.sources[0].shine(hitpoint = (-1,40,1))
    shining_c = beamLine.sources[0].shine(hitpoint = (1,40,-1))
    # Join 2 shining sources (or 1 shining in multiple directions)
    shining.concatenate(shining_b)
    shining.concatenate(shining_c)
    show = screen.expose(shining)
    plt.scatter(show.x, show.z)
    plt.show()


def get_photons():
    file = open('photons.csv', 'r')
    csv_reader = csv.reader(file, delimiter = '\t',\
                            quoting=csv.QUOTE_NONNUMERIC)
    photons = []
    for row in csv_reader:
        photons.append(row)
    return photons



def photo_positions(photons):

    # Prepare position containters
    xx = []
    zz = []
    for row in photons:
        xx.append(row[0])
        zz.append(row[1])

    print "Number of saved photons:", len(xx)
    plt.scatter(xx,zz)
    plt.show()

def photo_momenta(photons):

    # Prepare momenta containter
    xx = []
    zz = []
    aa = []
    bb = []
    cc = []
    for row in photons:
        xx.append(row[0])
        zz.append(row[1])
        aa.append(row[2])
        bb.append(row[3])
        cc.append(row[4])
    plt.quiver(xx, zz, aa, cc)
    plt.show()

if __name__ == '__main__':
    if False:
	photons = get_photons()
	photo_positions(photons)
	photo_momenta(photons)

    if True:
	view_source()
