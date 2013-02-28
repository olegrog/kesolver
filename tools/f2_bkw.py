#!/usr/bin/env python

import sys, math, array, struct, itertools
import json
import re

import numpy as np

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

if symmetry == "Cylindrical":
    v = float(data["gas"].get("v", "0."))
    dim = (2*rad, rad) 
    x = np.fromfunction(lambda i, j: v + (i + 0.5 - rad) * cut / rad,
                        (2 * rad, rad), dtype=float)
    y = np.fromfunction(lambda i, j: (j + 0.5) * cut / rad,
                        (2 * rad, rad), dtype=float)
    vol = y
    e = x*x + y*y
    
elif symmetry == "Cartesian":
    v = map(float, data["gas"].get("v", "( 0. 0. 0. )")[1:-1].split() )
    vx, vy, vz = v
    dim = (2*rad, 2*rad, 2*rad)
    x = np.fromfunction(lambda i, j, k: vx + (i + 0.5 - rad) * cut / rad,
                        (2 * rad, 2 * rad, 2 * rad), dtype=float)
    y = np.fromfunction(lambda i, j, k: vy + (j + 0.5 - rad) * cut / rad,
                        (2 * rad, 2 * rad, 2 * rad), dtype=float)
    z = np.fromfunction(lambda i, j, k: vz + (k + 0.5 - rad) * cut / rad,
                        (2 * rad, 2 * rad, 2 * rad), dtype=float)
    vol = np.ones_like(y)
    e = x*x + y*y + z*z

# calculate lambda

with open(maoutfilename, 'rb') as fd:
    size = int(fd.readline())
    theta, section = [], []
    for i in range(size):
        theta.append(float(fd.readline()))
    for i in range(size):
        section.append(float(fd.readline()))

#costheta = np.cos(theta)
#section = np.array(section)
#lamd = - math.pi / 2 * np.trapz(x=costheta, y=(section * (1 - sqr(costheta))), dx=1./size)

section = np.array(section)
lamd = np.trapz(x=theta, y=(section * cube(np.sin(theta))), dx=1./size)

with open(ffilename, 'rb') as fd:
    data = fd.read(4+4)
    i, size = struct.unpack('=ii', data)

    f = np.zeros_like(x)
    a = array.array('d')
    circl = e < cut*cut
    size = np.sum(circl)
    a.fromfile(fd, size)

    f[circl] = np.array(a)

# prepare function calculated by projection method

f = f / vol 
f = np.ma.masked_where(f == 0.0, f)

# prepare function calculated by theoretical formula

file_i = int(re.search(r"(\d+)[^/]*$", ffilename).groups()[0])
print file_i
time = file_i * timestep / math.sqrt(2) / math.pi

xi = 0.4
tau = 1. - xi * math.exp( - lamd * time )
e1 = 0.5 * e / tau

g = 1. / math.sqrt(cube(2 * math.pi * tau)) * np.exp(-e1) * ( 1 + (1 - tau) / tau * (e1 - 1.5) )

# prepare data for plotting

def p(symm, x, y, f):
    if symm == "Cartesian":
        xp = x[:, rad:, rad]
        yp = y[:, rad:, rad]
        fp = f[:, rad:, rad]
        return xp, yp, fp
    elif symm == "Cylindrical":
        return x, y, f 

# plot functions

surf = ax.plot_wireframe(*p(symmetry, x, y, f), color='b')
surg = ax.plot_wireframe(*p(symmetry, x, y, g), color='r')
#surg = ax.plot_wireframe(*p(symmetry, x, y, (f-g)), color='r')

with open("out1.txt", "w") as out:
    x1, y1, f1 = p(symmetry, x, y, f)
    for u, v in zip(x1[:, 0], f1[:, 0]):
        out.writelines(str(u) + ' ' + str(v) + '\n')

#with open("out2.txt", "w") as out:
#    x1, y1, f1 = p(symmetry, x, y, (f-g))
#    for u, v in zip(x1[:, 0], f1[:, 0]):
#        out.writelines(str(u) + ' ' + str(v) + '\n')

# estinate the error

err = math.sqrt(np.sum(sqr(f-g)) / np.sum(g))
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
