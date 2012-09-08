#!/usr/bin/env python

import sys, numpy
import re
import out2

def gmshTypeToVTKType(gmshtype):
    VTKTypes = [-1, -1, -1, -1, 10, 12, 13]
    return VTKTypes[gmshtype]

def listToStr(elm):
    return str(len(elm)) + ' ' + ' '.join( map(str, elm) )

def list3ToStr(l):
    return ' '.join( map(str, l) )

nodes, cells = out2.readNodesCells(sys.argv[1])
data = out2.readMacros(sys.argv[2], len(cells))

centernodes = []
for cell in cells:
    centernode = numpy.zeros(3)
    for vertex in cell.vertexes:
        centernode += nodes[vertex]
    centernodes.append(centernode / len(cell.vertexes))

with open(sys.argv[3], "w") as fd:
    fd.writelines("# vtk DataFile Version 2.0\n")
    fd.writelines("Velocity, MassFlux, HeatFlux, EnergyFlux.\n")    
    fd.writelines("ASCII\n")
    fd.writelines("DATASET UNSTRUCTURED_GRID\n")

    fd.writelines("POINTS %d float\n" % (len(nodes)+len(centernodes)) )
    for node in nodes:
        fd.writelines("%f %f %f\n" % (node[0], node[1], node[2]) )
    for centernode in centernodes:
        fd.writelines("%f %f %f\n" % (centernode[0], centernode[1], centernode[2]) )

    fd.writelines("CELLS %d %d\n" % ( len(cells), sum( [len(cell.vertexes)+1 for cell in cells] ) ))
    for cell in cells:
        fd.writelines(listToStr(cell.vertexes)+'\n')

    fd.writelines("CELL_TYPES %d\n" % len(cells))
    for cell in cells:
        fd.writelines(str(gmshTypeToVTKType(cell.type))+'\n')

    fd.writelines("CELL_DATA %d\n" % len(cells))

    for i, n in enumerate(data[0::8]):
        fd.writelines("SCALARS n%d float 1\n" % i)
        fd.writelines("LOOKUP_TABLE default\n")
        for density in n:
            fd.writelines(str(density)+'\n')

    for i, t in enumerate(data[4::8]):
        fd.writelines("SCALARS t%d float\n" % i)
        fd.writelines("LOOKUP_TABLE default\n")
        for temperature in t:
            fd.writelines(str(temperature)+'\n')

    fd.writelines("POINT_DATA %d\n" % (len(nodes)+len(centernodes)))
    for i, vx in enumerate(data[1::8]):
        v = data[1+i*8:4+i*8]
        fd.writelines("VECTORS v%d float\n" % i)
        for i in range(len(nodes)):
            fd.writelines('0. 0. 0.\n')
        for velocity in v:
            fd.writelines(list3ToStr(velocity)+'\n')



