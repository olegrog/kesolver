#!/usr/bin/env python

import sys, math, numpy, out2

O = numpy.array( [0., 0.0013, 0.001] )
u = numpy.array( [1., 0., 0.] )
v = numpy.array( [0., 0.997, 0.08] )

def intersect(O, u, v, cell, nodes):
	g = out2.gamma(nodes[cell.nodes[0]] - O, u, v)	
	for node in cell.nodes[1:]:
		if g * out2.gamma(nodes[node] - O, u, v) <= 0:
			return True
	return False

def angle(x, y):
	a = math.acos( numpy.dot(x, y) / 
			math.sqrt( numpy.dot(x, x) * numpy.dot(y, y) ) )
	if (numpy.cross(x, y) >= 0):
		return a
	else:
		return -a

def project(cell, nodes, O, u, v):
	ps = [numpy.array(p) for p in out2.cellCrossPlane(cell, nodes, O, u, v)]
	if ps != []:
		c = out2.center(ps)
		ps.sort( key = lambda x: angle( x-c, numpy.array( [1., 0.] ) ) )
		return ps


nodes, cells = out2.readNodesCells(sys.argv[1])
data = out2.readMacros(sys.argv[2], len(cells))

import matplotlib
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import pylab

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


width = 400
height = width/5
init_plot_params(width, height)

fig = pylab.figure()
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.07, right=0.99, bottom=0.15, top=0.95)
patches = []
colors = []

try:
    i = int(sys.argv[3])
except IndexError:
    i = 0

for cell, macro in zip(cells, numpy.array(data[i])):
#for cell, macro in zip(cells, numpy.array(data[9]) - numpy.array(data[17])):
#for cell, macro in zip(cells, numpy.array(data[i+8])/numpy.array(data[i])):
#for cell, macro in zip(cells, numpy.array(data[i])*numpy.array(data[i+4])):
	if intersect(O, u, v, cell, nodes):
		p = project(cell, nodes, O, u, v)		
		if p != None:
			xs = [x[0] for x in p] 
			ys = [x[1] for x in p] 

			polygon = Polygon(p, True)
			patches.append(polygon)

			colors.append(macro)


p = PatchCollection(patches, cmap = matplotlib.cm.jet, linewidths=1.)
p.set_edgecolor('face')
p.set_array(pylab.array(colors))
ax.add_collection(p)
ax.autoscale_view( tight=True )
pylab.colorbar(p)

pylab.show()

