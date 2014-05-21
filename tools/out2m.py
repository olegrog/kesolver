#!/usr/bin/env python

import sys, numpy, out2

nodes, cells = out2.readNodesCells(sys.argv[1])
volumes = [out2.cellVolume(cell, nodes) for cell in cells]
data    = out2.readMacros(sys.argv[2], len(cells))

indexes = [cell.phys_index for cell in cells]

print "sumvol = ", sum(volumes)

masses = {}
for index, volume, macro in zip(indexes, volumes, zip(*data)):
    mass = volume * numpy.array(macro[0::15])
    if index in masses:
        masses[index] += mass
    else:
        masses[index] = mass

with open(sys.argv[3], "w") as fd:
    for key, mass in zip(masses.keys(), masses.values()):
        if key >= 0:
            fd.writelines("%s %s\n" % ( key, ' '.join(map(str, mass)) ))

