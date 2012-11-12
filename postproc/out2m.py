#!/usr/bin/env python

import sys, numpy, out2

nodes, cells, indexes = out2.readNodesCellsIndexes(sys.argv[1])
volumes = [out2.cellVolume(cell, nodes) for cell in cells]
data = out2.readMacros(sys.argv[2], len(cells))

print "sumvol = ", sum(volumes)

filters = [lambda x: x[0] < -5.0, lambda x: x[0] > 5.0]
filters = []
if filters:
    def get_index(x, filters):
        for i, f in enumerate(filters):
            if f(x):
                return i
        return -1
    centers = [out2.center([nodes[i] for i in cell.vertexes]) for cell in cells]
    indexes = [get_index(x, filters) for x in centers]

masses = {}
for index, volume, macro in zip(indexes, volumes, zip(*data)):
    mass = volume * numpy.array(macro[0::8])
    if index in masses:
        masses[index] += mass
    else:
        masses[index] = mass

print "summass = ", sum(masses.values())

with open(sys.argv[3], "w") as fd:
    for key, mass in zip(masses.keys(), masses.values()):
        if key >= 0:
            fd.writelines("%d %s\n" % ( key, ' '.join(map(str, mass)) ))
