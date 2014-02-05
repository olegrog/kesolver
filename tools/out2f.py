#!/usr/bin/env python

import sys, math, array, struct, itertools
import json
import numpy as np

from mpl_toolkits.mplot3d.axes3d import Axes3D
import matplotlib.pyplot as plt

from kepy.ximesh import read_ximesh
from out2 import readNodesCells, center

def sqr(x):
    return x*x

def init_plot_params(columnwidth, height):
    fig_width_pt = columnwidth
    inches_per_pt = 1.0/72.27               # Convert pt to inch
    fig_width = fig_width_pt*inches_per_pt  # width in inches
    fig_height = height*inches_per_pt
    fig_size =  [fig_width,fig_height]
    params = {
                'backend': 'ps',
                'axes.labelsize': 10,
                'text.fontsize': 10,
                'legend.fontsize': 9,
                'xtick.labelsize': 7,
                'ytick.labelsize': 7,
                'text.usetex': False,
                'figure.figsize': fig_size,
                'font.family':'serif'
            }
    plt.rcParams.update(params)

width = 350
height = width * 0.75
init_plot_params(width, height)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
fig.subplots_adjust(left=-0.2, right=1.1, bottom=0, top=1.1)

# open .kei file
with open(sys.argv[1], 'rb') as fd:
    data = json.load(fd)

nodes, cells = readNodesCells(sys.argv[1])

symmetry, rad, circl, xyz, vol, r, d3v = read_ximesh(data)

with open(sys.argv[2], 'rb') as fd:
    while fd:
        data = fd.read(4+4)
        i, size = struct.unpack('=ii', data)

        print i, "size (file) = ", size

        print center([nodes[n] for n in cells[i].nodes])


        f = np.zeros_like(r)
        f.fill(-1e-10)
        a = array.array('d')
        size = np.sum(circl)
        print "size (ximesh) = ", size
        a.fromfile(fd, size)

        if i == int(sys.argv[3]):
            f[circl] = np.array(a)

            f = f / vol

            if symmetry == "Cartesian":
                xp = xyz[0][:, rad:, rad]
                yp = xyz[1][:, rad:, rad]
                fp = f[:, rad:, rad]
            elif symmetry == "Cylindrical":
                xp, yp, fp = xyz[0], xyz[1], f

            print fp < 0.

            fp = np.ma.masked_where(fp < 0., fp)
            surf = ax.plot_wireframe(yp, xp, fp, color='k',linewidth=0.5, rstride=2, cstride=2)

            break

ax.set_xlabel(r'$\xi_y$')
ax.set_ylabel(r'$\xi_x$')

for axis in ax.w_xaxis, ax.w_yaxis, ax.w_zaxis:
    axis.pane.set_visible(False)
    axis.gridlines.set_visible(False)
    axis.set_rotate_label(False)

for elt in ax.w_zaxis.get_ticklines() + ax.w_zaxis.get_ticklabels():
    elt.set_visible(False)

ax.view_init(20, -5)
fig.savefig("f3d.eps")

plt.show()

