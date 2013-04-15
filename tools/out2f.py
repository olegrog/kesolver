#!/usr/bin/env python

import sys, math, array, struct, itertools
import json
import numpy as np

from mpl_toolkits.mplot3d.axes3d import Axes3D
import matplotlib.pyplot as plt

from kepy.ximesh import read_ximesh

def sqr(x):
    return x*x

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# open .kei file
with open(sys.argv[1], 'rb') as fd:
    data = json.load(fd)

symmetry, rad, circl, xyz, vol, r, d3v = read_ximesh(data)

with open(sys.argv[2], 'rb') as fd:
    while fd:
        data = fd.read(4+4)
        i, size = struct.unpack('=ii', data)

        print i, "size (file) = ", size

        f = np.zeros_like(r)
        a = array.array('d')
        size = np.sum(circl)
        print "size (ximesh) = ", size
        a.fromfile(fd, size)

        if i == int(sys.argv[3]):
            f[circl] = np.array(a)

            f = f / vol
            f = np.ma.masked_where(f == 0.0, f)

            if symmetry == "Cartesian":
                xp = xyz[0][:, rad:, rad]
                yp = xyz[1][:, rad:, rad]
                fp = f[:, rad:, rad]
            elif symmetry == "Cylindrical":
                xp, yp, fp = xyz[0], xyz[1], f
            surf = ax.plot_wireframe(xp, yp, fp, color='k')

            break

ax.set_xlabel(r'$\xi_x$')
ax.set_ylabel(r'$\xi_y$')

for axis in ax.w_xaxis, ax.w_yaxis, ax.w_zaxis:
    axis.pane.set_visible(False)
    axis.gridlines.set_visible(False)
    axis.set_rotate_label(False)

for elt in ax.w_zaxis.get_ticklines() + ax.w_zaxis.get_ticklabels():
    elt.set_visible(False)

plt.show()
