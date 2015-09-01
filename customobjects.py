import numpy as np
import xrt.backends.raycing.oes as roe
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches

class CustomShape(roe.OE):
    def __init__(self, *args, **kwargs):
        self.R = kwargs.pop('R', 0.5)
        roe.OE.__init__(self, *args, **kwargs)

    def local_n(self, x, y):
        # x
        a = 0.0
        # y
        b = 4.0
        # z
        c = 1.0
        norm = (a**2 + c**2)**0.5
        return [a/norm, b/norm, c/norm]

    def local_z(self, x, y):
        return -2. + 0.2*(y - self.limPhysY[0])

if __name__ == '__main__':
    verts = [
        (0., 0.), # left, bottom
        (0., 3.), # left, top
        (1., 3.), # right, top
        (1., 1.),
        (2., 1.),
        (2., 3.),
        (3., 3.),
        (3., 0.),
        (2., 0.), # right, bottom
        (2., -2),
        (1., -2),
        (1., 0.),
        (0., 0.), # ignored
        ]

    verts = [(v[0] - 1.5, v[1]) for v in verts]
    verts = [(0.1*v[0], 0.1*v[1]) for v in verts]

    path = Path(verts)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    patch = patches.PathPatch(path, facecolor='orange', lw=2)
    ax.add_patch(patch)
    ax.set_xlim(-4,4)
    ax.set_ylim(-4,4)
    plt.show()
