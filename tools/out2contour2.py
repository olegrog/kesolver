#!/usr/bin/env python

import sys, math, numpy, out2

O = numpy.array( [0.0012, 0.0013, 0.001] )
u = numpy.array( [1., 0., 0.] )
v = numpy.array( [0., 1.0, 0.0] )

#O = numpy.array( [0., 0.0013, 0.001] )
#u = numpy.array( [1., 0., 0.] )
#v = numpy.array( [0., 0.8, 0.6] )

nodes, cells = out2.readNodesCells(sys.argv[1])
data = out2.readMacros(sys.argv[2], len(cells))

verges = out2.cellsToVerges(cells)
points = out2.intersectEdges(verges, cells, nodes, O, u, v)
colors = numpy.zeros(len(points))
strikes = numpy.zeros(len(points))

import matplotlib
from matplotlib.ticker import MaxNLocator
import pylab
import matplotlib.tri as tri

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
                'xtick.labelsize': 8,
                'ytick.labelsize': 8,
                'text.usetex': False,
                'figure.figsize': fig_size,
                'font.family':'serif'
            }
    pylab.rcParams.update(params)

def make_cmap():
    cdict = {'red'   : ((0.000, 0.00, 0.00),
                        (0.125, 0.15, 0.15),
                        (0.250, 0.30, 0.30),
                        (0.375, 0.60, 0.60),
                        (0.500, 1.00, 1.00),
                        (0.625, 0.90, 0.90),
                        (0.750, 0.90, 0.90),
                        (0.875, 0.90, 0.90),
                        (1.000, 1.00, 1.00)),
             
             'green' : ((0.000, 0.00, 0.00),
                        (0.125, 0.15, 0.15),
                        (0.250, 0.15, 0.15),
                        (0.375, 0.20, 0.20),
                        (0.500, 0.25, 0.25),
                        (0.625, 0.50, 0.50),
                        (0.750, 0.75, 0.75),
                        (0.875, 0.90, 0.90),
                        (1.000, 1.00, 1.00)),
             
             'blue'  : ((0.000, 0.00, 0.00),
                        (0.125, 0.50, 0.50),
                        (0.250, 0.75, 0.75),
                        (0.375, 0.50, 0.50),
                        (0.500, 0.15, 0.15),
                        (0.625, 0.00, 0.00),
                        (0.750, 0.10, 0.10),
                        (0.875, 0.50, 0.50),
                        (1.000, 1.00, 1.00))
             }
    cm = matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 1024, gamma=.5)
    return cm

width = 225 
height = width*3/6
init_plot_params(width, height)

fig = pylab.figure()
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.06, right=0.94, bottom=0.05, top=0.95)

try:
    i = int(sys.argv[3])
except IndexError:
    i = 0

x, y = zip(*points)
x, y = numpy.array(x), numpy.array(y)
#trs = tri.Triangulation(x, y)
trs = []

for cell, macro in zip(cells, numpy.array(data[i])):
#for cell, macro in zip(cells, numpy.array(data[i]) * numpy.array(data[i+4])):
#for cell, macro in zip(cells, numpy.array(data[i+8])/numpy.array(data[i])):
#for cell, macro in zip(cells, numpy.array(data[i])*numpy.array(data[i+4])/numpy.array(data[i+8])/numpy.array(data[i+12])):
    if cell.points:
        trs += out2.triangles(cell.points, points)
        for p in cell.points:
            colors[p]  += macro
            strikes[p] += 1

colors /= strikes

#size = len(x)
#x = numpy.concatenate( (x, x, x) )
#y = numpy.concatenate( (y, -y, 2-y) )
#colors = numpy.concatenate( (colors, colors, colors) )

#newtrs = [[k + size for k in tr] for tr in trs]
#newtrs2 = [[k + 2*size for k in tr] for tr in trs]
#trs += newtrs+newtrs2

cm = make_cmap()

ax.autoscale_view( tight=True )
#co = ax.tricontourf(-numpy.array(x), 1 - numpy.array(y), trs, colors, 20)
#co = ax.tricontourf(13-numpy.array(x), numpy.array(y), trs, colors, 20)
co = ax.tricontourf(numpy.array(x), numpy.array(y), trs, colors, 20, cmap = cm)

bar = pylab.colorbar(co, orientation='horizontal')

xmajorLocator   = MaxNLocator(10)
ax.xaxis.set_major_locator(xmajorLocator)

ymajorLocator   = MaxNLocator(5)
ax.yaxis.set_major_locator(ymajorLocator)

#barmajorLocator   = MaxNLocator(5)
#bar.ax.xaxis.set_major_locator(barmajorLocator)

pylab.savefig("field.eps")
pylab.show()

