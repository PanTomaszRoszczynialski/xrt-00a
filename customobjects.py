import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches

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
verts = [(0.3*v[0], 0.3*v[1]) for v in verts]

path = Path(verts)

fig = plt.figure()
ax = fig.add_subplot(111)
patch = patches.PathPatch(path, facecolor='orange', lw=2)
ax.add_patch(patch)
ax.set_xlim(-4,4)
ax.set_ylim(-4,4)
plt.show()
