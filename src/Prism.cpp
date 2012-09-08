#include <iostream>
#include <fstream>

#include "Prism.hpp"
#include "Tetrahedron.hpp"

void Prism::calculateVolume()
{
	V = tetrahedronVolume(vertex[0], vertex[1], vertex[2], vertex[3]) + 
		tetrahedronVolume(vertex[1], vertex[3], vertex[4], vertex[5]) + 
		tetrahedronVolume(vertex[1], vertex[2], vertex[3], vertex[5]);
}


void Prism::calculateCenter()
{
	double V1 =	tetrahedronVolume(vertex[0], vertex[1], vertex[2], vertex[3]);
	V3d c1 =	tetrahedronCenter(vertex[0], vertex[1], vertex[2], vertex[3]);

	double V2 =	tetrahedronVolume(vertex[1], vertex[3], vertex[4], vertex[5]);
	V3d c2 =	tetrahedronCenter(vertex[1], vertex[3], vertex[4], vertex[5]);

	double V3 =	tetrahedronVolume(vertex[1], vertex[2], vertex[3], vertex[5]);
	V3d c3 =	tetrahedronCenter(vertex[1], vertex[2], vertex[3], vertex[5]);

	center = (c1 * V1 + c2 * V2 + c3 * V3) / (V1 + V2 + V3);
}


Prism::Prism()
{
	numberOfVertex = 6;
	numberOfEdges = 9;
	numberOfFacets = 5;
}
