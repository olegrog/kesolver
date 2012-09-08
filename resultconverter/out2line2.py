#!/usr/bin/env python

import sys, math, numpy
import pylab
import out2
from math import pi, sin, cos

nodes, cells = out2.readNodesCells(sys.argv[1])
data = out2.readMacros(sys.argv[2], len(cells))

def make_snake(m, n, l, R):
    xs = []
    ys = []
    for i in range(m):
        x0 = i * 4 * R
        xs += [x0                for j in range(n)]
        ys += [(j + 0.5) / n * l for j in range(n)]

        x0 = i * 4 * R + R
        xs += [x0 - R * cos((j + 0.5) / n * pi) for j in range(n)]
        ys += [l  + R * sin((j + 0.5) / n * pi) for j in range(n)]
    
        x0 = i * 4 * R + 2 * R
        xs += [x0                      for j in range(n)]
        ys += [(1 - (j + 0.5) / n) * l for j in range(n)]

        x0 = i * 4 * R + 3 * R
        xs += [x0 - R * cos((j + 0.5) / n * pi) for j in range(n)]
        ys += [   - R * sin((j + 0.5) / n * pi) for j in range(n)]

    return xs, ys
    
xs, ys = make_snake(5, 10, 4, 1.3)
path = [numpy.array([x+0.001, y+0.001, -0.001]) for x, y in zip(xs, ys)]

nxs, nys, nzs = zip(*nodes)
xmin, xmax = min(nxs), max(nxs)
ymin, ymax = min(nys), max(nys)
zmin, zmax = min(nzs), max(nzs)
print zmin, zmax

nx = 100
ny = 20
nz = 5 
grid_x, grid_y, grid_z = numpy.mgrid[xmin:xmax:(nx*1j), ymin:ymax:(ny*1j), zmin:zmax:(nz*1j)]
m = [[[[] for k in range(nz)] for j in range(ny)] for i in range(nx)]
for i, cell in enumerate(cells):
    xs, ys, zs = zip(*[nodes[v] for v in cell.vertexes])
    x1, x2 = min(xs), max(xs)
    y1, y2 = min(ys), max(ys)
    z1, z2 = min(zs), max(zs)
    i1 = max(0,  int( (x1 - xmin) / (xmax - xmin) * (nx - 1) ))
    i2 = min(nx, int( (x2 - xmin) / (xmax - xmin) * (nx - 1) ) + 2)
    j1 = max(0,  int( (y1 - ymin) / (ymax - ymin) * (ny - 1) ))
    j2 = min(ny, int( (y2 - ymin) / (ymax - ymin) * (ny - 1) ) + 2)
    k1 = max(0,  int( (z1 - zmin) / (zmax - zmin) * (nz - 1) ))
    k2 = min(nz, int( (z2 - zmin) / (zmax - zmin) * (nz - 1) ) + 2)
    ii = [ (xi, yi, zi) for xi in range(i1, i2) for yi in range(j1, j2) for zi in range(k1, k2)]
    for xi, yi, zi in ii:
        m[xi][yi][zi].append(i)

size = len(cells)
celli = []
ti    = []
for i, p in enumerate(path):
    x, y, z = p
    xi = min(nx - 1, max(0, int( (x - xmin) / (xmax - xmin) * (nx - 1) + 0.5)))
    yi = min(ny - 1, max(0, int( (y - ymin) / (ymax - ymin) * (ny - 1) + 0.5)))
    zi = min(nz - 1, max(0, int( (z - zmin) / (zmax - zmin) * (nz - 1) + 0.5)))
    for j in m[xi][yi][zi]:
        if out2.inCell(cells[j], nodes, p):
            celli.append(j)
            ti.append(i)
            break
print celli

at = lambda array, indexes : [array[i] for i in indexes] 
atdata = [at(macro, celli) for macro in data]

toprint = []
for i, macro in enumerate(atdata):
    pylab.plot(ti, macro, label=str(i))
    toprint.append(macro)

pylab.grid(True)
pylab.show()

if True:
    toprint = zip(*toprint)
    with open("line.txt", "w") as fd:
        for x, y in zip(ti, toprint):
            fd.writelines("%f %s\n" % (
                x, 
                ' '.join(map(str, y))
            ))
