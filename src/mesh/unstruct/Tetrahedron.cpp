#include <iostream>
#include <fstream>

#include "Tetrahedron.hpp"

void Tetrahedron::calculateVolume()
{
	V = tetrahedronVolume(vertex[0], vertex[1], vertex[2], vertex[3]);
}


void Tetrahedron::calculateCenter()
{
	center = tetrahedronCenter(vertex[0], vertex[1], vertex[2], vertex[3]);
}


Tetrahedron::Tetrahedron()
{
	numberOfVertex = 4;
	numberOfEdges  = 6;
	numberOfFacets = 4;
}
