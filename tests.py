import csv
import numpy as np
import matplotlib.pyplot as plt

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
	if len(row) == 5:
	    xx.append(row[0])
	    zz.append(row[1])
	else:
	    print 'Bad row length, check writing mechanism'

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
	if len(row) == 5:
	    xx.append(row[0])
	    zz.append(row[1])
	    aa.append(row[2])
	    bb.append(row[3])
	    cc.append(row[4])
	else:
	    print 'Bad row length, check writing mechanism'

    plt.quiver(xx, zz, aa, cc)
    plt.show()

if __name__ == '__main__':
    photons = get_photons()
    photo_positions(photons)
    photo_momenta(photons)
