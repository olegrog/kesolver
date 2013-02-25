#!/usr/bin/env python

import sys, math, numpy, array, struct, itertools
import json
import re

from mpl_toolkits.mplot3d.axes3d import Axes3D
import matplotlib.pyplot as plt

from scipy.integrate import quad

def sqr(x):
    return x * x

def cube(x):
    return x * x * x

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

keifilename   = sys.argv[1]
maoutfilename = sys.argv[2]
ffilename     = sys.argv[3]
timestep      = float(sys.argv[4])

# open .kei file
with open(keifilename, 'rb') as fd:
    data = json.load(fd)

cut = float(data["gas"]["cut"])
rad = int(data["gas"]["rad"])
symmetry = data["gas"]["symmetry"]

print rad

dim = (2*rad, rad) 
x = numpy.fromfunction(lambda i, j: (i + 0.5 - rad) * cut / rad, dim, dtype=float)
r = numpy.fromfunction(lambda i, j: (j + 0.5) * cut / rad, dim, dtype=float)

# calculate lambda

with open(maoutfilename, 'rb') as fd:
    size = int(fd.readline())
    theta, section = [], []
    for i in range(size):
        theta.append(float(fd.readline()))
    for i in range(size):
        section.append(float(fd.readline()))

#costheta = numpy.cos(theta)
#section = numpy.array(section)
#lamd = - math.pi / 2 * numpy.trapz(x=costheta, y=(section * (1 - sqr(costheta))), dx=1./size)

section = numpy.array(section)
lamd = numpy.trapz(x=theta, y=(section * cube(numpy.sin(theta))), dx=1./size)

print lamd


with open(ffilename, 'rb') as fd:
    data = fd.read(4+4)
    i, size = struct.unpack('=ii', data)

    f = numpy.zeros_like(x)
    a = array.array('d')
    circl = x*x + r*r < cut*cut
    size = numpy.sum(circl)
    a.fromfile(fd, size)

    f[circl] = numpy.array(a)

# plot function calculated by projection method

f = f / r
f = numpy.ma.masked_where(f == 0.0, f)
surf = ax.plot_wireframe(x, r, f, color='b')

# plot function calculated by theoretical formula

file_i = int(re.search(r"(\d+)[^/]*$", ffilename).groups()[0])
print file_i
time = file_i * timestep / math.sqrt(2) / math.pi

xi = 0.4
tau = 1. - xi * math.exp( - lamd * time )
e = 0.5 * ( sqr(x) + sqr(r) ) / tau

g = 1. / math.sqrt(cube(2 * math.pi * tau)) * numpy.exp(-e) * ( 1 + (1 - tau) / tau * (e - 1.5) )

surg = ax.plot_wireframe(x, r, g, color='r')

#surg = ax.plot_wireframe(x, r, (f-g), color='r')

err = math.sqrt(numpy.sum(sqr(f-g)) / numpy.sum(g))
print err

ax.set_xlabel(r'$\xi_y$')
ax.set_ylabel(r'$\xi_x$')

for axis in ax.w_xaxis, ax.w_yaxis, ax.w_zaxis:
    axis.pane.set_visible(False)
    axis.gridlines.set_visible(False)
    axis.set_rotate_label(False)

#for elt in ax.w_zaxis.get_ticklines() + ax.w_zaxis.get_ticklabels():
#    elt.set_visible(False)

plt.show()
