#!/usr/bin/env python

import sys, math, numpy, os
import pylab

from out2 import *

O = numpy.array( [0.012, 0.011, 0.01] )
u = numpy.array( [1.,  0.,    0.] )

show = True
save = True

nodes, cells = readNodesCells(sys.argv[1])

points = []
inds = []

print len(cells)

for i, cell in enumerate(cells):
    point = intersect(O, u, cell, nodes)
    if point:
        points.append(point)
        inds.append(i)

points, inds = zip(*sorted(zip(points, inds), key = lambda pair: pair[0]))

at = lambda array, indexes : [array[i] for i in indexes] 

fig = pylab.figure()

reads = {'.txt': readMacros, '.bin': readF}

for filename in sys.argv[2:]:

    toprint = []

    data = readMacros(filename, len(cells))
        
    atdata = [at(macro, inds) for macro in data]

    for i, macro in enumerate(atdata):
        pylab.plot(points, macro, label=str(i))
        toprint.append(macro)

    if save:
        toprint = zip(*toprint)
        outfilename = os.path.splitext(filename)[0] + '.line'
        print outfilename
        with open(outfilename, "w") as fd:
            for x, y in zip(points, toprint):
                fd.writelines("%f %s\n" % (
                    x, 
                    ' '.join(map(str, y))
                ))

if show:
#   pylab.legend(frameon = False)
    pylab.legend()
    pylab.grid(True)
    pylab.show()

