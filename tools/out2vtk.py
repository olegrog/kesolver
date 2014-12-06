#!/usr/bin/env python

import sys, numpy
import re
import out2
import numpy as np

N = 15 # number of macroscopic variables

def gmshTypeToVTKType(gmshtype):
    VTKTypes = [-1, -1, -1, -1, 10, 12, 13]
    return VTKTypes[gmshtype]

def listToStr(elm):
    return str(len(elm)) + ' ' + ' '.join( map(str, elm) )

def list3ToStr(l):
    return ' '.join( map(str, l) )

def toVec(data, pos):
    return map(list, zip(*data[pos : pos + 3]))

def writeScalar(name, data, pos):
    for i, scalars in enumerate(data[pos::N]):
        fd.writelines("SCALARS " + name + "%d float\n" % i)
        fd.writelines("LOOKUP_TABLE default\n")
        for scalar in scalars:
            fd.writelines(str(scalar) + '\n')

def writeVector(name, data, pos):
    for i, vectors_x in enumerate(data[pos::N]):
        vectors = toVec(data, pos+i*N)
        fd.writelines("VECTORS " + name + "%d float\n" % i)
        for vector in vectors:
            fd.writelines(list3ToStr(vector) + '\n')

nodes, cells = out2.readNodesCells(sys.argv[1])
data = out2.readMacros(sys.argv[2], len(cells))

centernodes = []
for cell in cells:
    centernode = numpy.zeros(3)
    for vertex in cell.nodes:
        centernode += nodes[vertex]
    centernodes.append(centernode / len(cell.nodes))

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

    fd.writelines("CELLS %d %d\n" % ( len(cells), sum( [len(cell.nodes)+1 for cell in cells] ) ))
    for cell in cells:
        fd.writelines(listToStr(cell.nodes) + '\n')

    fd.writelines("CELL_TYPES %d\n" % len(cells))
    for cell in cells:
        fd.writelines(str(gmshTypeToVTKType(cell.type)) + '\n')

    fd.writelines("CELL_DATA %d\n" % len(cells))

    writeScalar('Density', data, 0)
    writeVector('Velocity', data, 1)
    writeScalar('Temperature', data, 4)
    writeVector('Temperatures', data, 5)
    writeVector('HeatFlux', data, 8)
    writeVector('ShearStress', data, 11)
    writeScalar('H-function', data, 14)
