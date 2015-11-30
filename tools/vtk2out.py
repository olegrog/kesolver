#!/usr/bin/env python

import numpy as np
import sys, os, vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

vtkfile = sys.argv[1]
mfile = sys.argv[2]
total_rho = float(sys.argv[3]) if len(sys.argv) > 3 else 1

def get_cell_data(out, name):
    vtk_array = out.GetCellData().GetArray(name)
    return vtk_to_numpy(vtk_array)

vec2str = lambda l: '( ' + ' '.join( map(str, l) ) + ' )'

### Read a VTK-file
reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName(vtkfile)
reader.Update()
out = reader.GetOutput()
N = out.GetNumberOfCells()

### Read fields
rho = get_cell_data(out, 'rho')
T = get_cell_data(out, 'T')
U = get_cell_data(out, 'U')
vol = np.array([ vtk.vtkMeshQuality.HexVolume(out.GetCell(i)) for i in xrange(N) ])

### Correct rho
print 'Initial rho = %.7g' % (np.sum(rho*vol) / np.sum(vol))
rho *= total_rho / np.sum(rho*vol) * np.sum(vol)
print 'Resulting rho = %.7g' % (np.sum(rho*vol) / np.sum(vol))

zs, zv = 0, vec2str([0,0,0])

# NB: np.sqrt(2) Sone --> Tcheremissine
with open(mfile, "w") as f:
    for i in xrange(len(T)):
        f.write("%d [ %g %s %g %s %s %s %g ]\n" % \
            (i, rho[i], vec2str(U[i]*np.sqrt(2)), T[i], zv, zv, zv, zs))  
    
