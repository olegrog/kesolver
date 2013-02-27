#!/usr/bin/env python

import sys, math, numpy, array, struct, itertools
import json

from mpl_toolkits.mplot3d.axes3d import Axes3D
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# open .kei file
with open(sys.argv[1], 'rb') as fd:
    data = json.load(fd)

cut = float(data["gas"]["cut"])
rad = int(data["gas"]["rad"])
symmetry = data["gas"]["symmetry"]

print rad

if symmetry == "Cylindrical":
    dim = (2*rad, rad) 
    x = numpy.fromfunction(lambda i, j: i - rad + 0.5, dim)
    r = numpy.fromfunction(lambda i, j: j       + 0.5, dim)
    vol = r
    e = x*x + r*r
    
elif symmetry == "Cartesian":
    dim = (2*rad, 2*rad, 2*rad)
    x = numpy.fromfunction(lambda i, j, k: i - rad + 0.5, dim)
    y = numpy.fromfunction(lambda i, j, k: j - rad + 0.5, dim)
    z = numpy.fromfunction(lambda i, j, k: k - rad + 0.5, dim)
    vol = numpy.ones_like(y)
    e = x*x + y*y + z*z


with open(sys.argv[2], 'rb') as fd:
    while fd:
        data = fd.read(4+4)
        i, size = struct.unpack('=ii', data)

        print i, size

        f = numpy.zeros_like(x)
        a = array.array('d')
        circl = e < rad*rad
        size = numpy.sum(circl)
        a.fromfile(fd, size)

        if i == int(sys.argv[3]):
            f[circl] = numpy.array(a)

            f = f / vol
            f = numpy.ma.masked_where(f == 0.0, f)

            if symmetry == "Cartesian":
                xp = x[:, rad:, rad]
                yp = y[:, rad:, rad]
                fp = f[:, rad:, rad]
                surf = ax.plot_wireframe(xp, yp, fp, color='k')
            elif symmetry == "Cylindrical":
                surf = ax.plot_wireframe(x, r, f, color='k')

            break

ax.set_xlabel(r'$\xi_y$')
ax.set_ylabel(r'$\xi_x$')

for axis in ax.w_xaxis, ax.w_yaxis, ax.w_zaxis:
    axis.pane.set_visible(False)
    axis.gridlines.set_visible(False)
    axis.set_rotate_label(False)

for elt in ax.w_zaxis.get_ticklines() + ax.w_zaxis.get_ticklabels():
    elt.set_visible(False)

plt.show()
