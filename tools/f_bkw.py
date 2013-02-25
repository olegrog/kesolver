#!/usr/bin/env python

import sys, math, numpy, array, struct, itertools
import json

from mpl_toolkits.mplot3d.axes3d import Axes3D
import matplotlib.pyplot as plt

def sqr(x):
    return x * x

def cube(x):
    return x * x * x

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# open .kei file
with open(sys.argv[1], 'rb') as fd:
    data = json.load(fd)

cut = float(data["gas"]["cut"])
rad = int(data["gas"]["rad"])
symmetry = data["gas"]["symmetry"]

print rad

dim = (2*rad, rad) 
x = numpy.fromfunction(lambda i, j: i - rad + 0.5, dim)
r = numpy.fromfunction(lambda i, j: j       + 0.5, dim)

fs = []

for filename in sys.argv[2:4]:
    with open(filename, 'rb') as fd:
        while fd:
            data = fd.read(4+4)
            i, size = struct.unpack('=ii', data)

            print i, size

            f = numpy.zeros_like(x)
            a = array.array('d')
            circl = x*x + r*r < rad*rad
            size = numpy.sum(circl)
            a.fromfile(fd, size)

            if i == int(sys.argv[4]):
                f[circl] = numpy.array(a)

                f = f / r
                fs.append(f)

                break

f = fs[1] - fs[0]

f = numpy.ma.masked_where(f == 0.0, f)
surf = ax.plot_wireframe(x, r, f, color='b')

xi = 0.4
tau = 1. - xi
e = 0.5 * ( sqr(x) + sqr(r) ) / tau
g = 1. / 15 / tau / math.sqrt(cube(2 * math.pi * tau)) * (1. / tau - 1) * numpy.exp(-e) * (sqr(e-1.5) + 1.5 - e)

surg = ax.plot_wireframe(x, r, g, color='r')

ax.set_xlabel(r'$\xi_y$')
ax.set_ylabel(r'$\xi_x$')

for axis in ax.w_xaxis, ax.w_yaxis, ax.w_zaxis:
    axis.pane.set_visible(False)
    axis.gridlines.set_visible(False)
    axis.set_rotate_label(False)

for elt in ax.w_zaxis.get_ticklines() + ax.w_zaxis.get_ticklabels():
    elt.set_visible(False)

plt.show()
