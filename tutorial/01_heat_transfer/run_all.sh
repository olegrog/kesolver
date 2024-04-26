#!/bin/bash

# Remove comments from a JSON file to generate a KESolver problem file
sed 's_ *//.*__' heat.kep_with_comments > heat.kep

# Generate a mesh using GMSH producing file `heat.msh`
# Note that we use a specific mesh format version 2.2
printf "Generating mesh..."
gmsh -3 -format msh22 heat.geo > log.gmsh
echo "done."

# Convert GMSH format to KESolver mesh format
../../tools/msh2kem.py heat.msh heat.kem

# Generate an input KESolver file from the problem and mesh files
../../tools/attach_mesh.py heat.kep heat.kem heat.kei

# Run solver
if [[ -d result ]]; then
    echo "Solver has been already run!"
else
    printf "Running solver..."
    ../../src/kes heat.kei > log.kes 2>&1
    echo "done."
fi

# Convert all files to the VTK format
if [[ -d vtk ]]; then
    echo "Result are already converted!"
else
    mkdir vtk
    printf "Converting results to the VTK format..."
    for out in $(ls result | grep '.txt'); do
        ../../tools/out2vtk.py heat.kei result/$out vtk/${out%txt}vtk
    done
    echo "done."
fi
