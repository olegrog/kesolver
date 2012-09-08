#!/usr/bin/env python

import sys, math, numpy, out2

O = numpy.array( [0., 0.00013, 0.0001] )
u = numpy.array( [1., 0., 0.] )
v = numpy.array( [0., 0.997, 0.08] )

#O = numpy.array( [0., 0., -0.001] )
#u = numpy.array( [1., 0., 0.] )
#v = numpy.array( [0., 1., 0.] )

nodes, cells  = out2.readNodesCells(sys.argv[1])
data = out2.readMacros(sys.argv[2], len(cells))

verges = out2.cellsToVerges(cells)
points = out2.intersectEdges(verges, cells, nodes, O, u, v)

import matplotlib
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
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


width = 225
height = width/2
init_plot_params(width, height)

fig = pylab.figure()

try:
    i = int(sys.argv[3])
except IndexError:
    i = 0

x, y = zip(*points)
x, y = numpy.array(x), numpy.array(y)
#trs = tri.Triangulation(x, y)
trs = []
colors  = numpy.array([numpy.zeros(3) for j in range(len(points))])
strikes = numpy.zeros(len(points))

for cell, line in zip(cells, zip(*data)):
    w = numpy.array(line[i:i+3])
    #w = numpy.array(line[17:20]) - numpy.array(line[9:12])
    #w = numpy.array(line[9:12]) - numpy.array(line[1:4])
    #w = numpy.array(line[17:20]) - numpy.array(line[1:4])
    if cell.points:
        trs += out2.triangles(cell.points, points)
        for p in cell.points:
            colors[p]  += w
            strikes[p] += 1

colors = [out2.projToPlane(clr / str, u, v) for clr, str in zip(colors, strikes)]
u, v = zip(*colors)
u, v = numpy.array(u), numpy.array(v)

"""
size = len(x)
x = numpy.concatenate( (x, x) )
y = numpy.concatenate( (y, -y) )
points = [numpy.array( (x1, y1) ) for x1, y1 in zip(x, y)]
u = numpy.concatenate( (u, u) )
v = numpy.concatenate( (v, -v) )

newtrs = [[k + size for k in tr] for tr in trs]
trs += newtrs
"""

z = numpy.sqrt(u*u + v*v)

ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.06, right=0.92, bottom=0.12, top=0.95)
cont = ax.tricontourf(x, y, trs, z, 25, cmap=pylab.cm.jet)
bar = pylab.colorbar(cont, orientation='horizontal')

xmin, xmax = numpy.min(x), numpy.max(x)
ymin, ymax = numpy.min(y), numpy.max(y)

print xmin, xmax
print ymin, ymax

nx = 300
ny = 100
grid_x, grid_y = numpy.mgrid[xmin:xmax:(nx*1j), ymin:ymax:(ny*1j)]
spx = numpy.linspace(xmin, xmax, nx)
spy = numpy.linspace(ymin, ymax, ny)

m  = numpy.ones_like(grid_x)
for tr in trs:
    xs = x[tr]
    ys = y[tr]
    x1, x2 = min(xs), max(xs)
    y1, y2 = min(ys), max(ys)
    i1 = max(0,      int( (x1 - xmin) / (xmax - xmin) * (nx - 1) ))
    i2 = min(nx - 1, int( (x2 - xmin) / (xmax - xmin) * (nx - 1) ) + 2)
    j1 = max(0,      int( (y1 - ymin) / (ymax - ymin) * (ny - 1) ))
    j2 = min(ny - 1, int( (y2 - ymin) / (ymax - ymin) * (ny - 1) ) + 2)
    ii = [ (xi, yi) for xi in range(i1, i2) for yi in range(j1, j2) ]
    for xi, yi in ii:
        if out2.inTr(numpy.array( (spx[xi], spy[yi]) ), tr, points):
            m[xi][yi] = 0

from scipy.interpolate import griddata
from matplotlib.ticker import MaxNLocator

grid_u = griddata(points, u, (grid_x, grid_y), method='linear')
grid_v = griddata(points, v, (grid_x, grid_y), method='linear')

grid_uv = numpy.sqrt(grid_u*grid_u + grid_v*grid_v)
max_uv = numpy.nanmax(grid_uv)
print max_uv
w = 10.;
m = numpy.logical_or(m, grid_uv < max_uv / w)
grid_u = numpy.ma.array(grid_u, mask=m)

grid_u  = grid_u.transpose()
grid_v  = grid_v.transpose()
grid_uv = grid_uv.transpose()

ax.streamplot(spx, spy, grid_u, grid_v, color='k', density=1,
             linewidth = 0.5, arrowstyle='-|>', arrowsize=0.5)
#ax.contourf(grid_x, grid_y, m)

ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)

xmajorLocator   = MaxNLocator(7)
ax.xaxis.set_major_locator(xmajorLocator)

ymajorLocator   = MaxNLocator(5)
ax.yaxis.set_major_locator(ymajorLocator)

barmajorLocator   = MaxNLocator(5)
bar.ax.xaxis.set_major_locator(barmajorLocator)

pylab.savefig("stream.eps")
pylab.show()


