#include <iostream>
#include <fstream>

#include "Hexahedron.hpp"
#include "Tetrahedron.hpp"

void Hexahedron::calculateVolume()
{
	V = tetrahedronVolume(vertex[0], vertex[1], vertex[3], vertex[4]) + 
		tetrahedronVolume(vertex[1], vertex[2], vertex[3], vertex[6]) + 
		tetrahedronVolume(vertex[1], vertex[4], vertex[5], vertex[6]) +  
		tetrahedronVolume(vertex[3], vertex[4], vertex[6], vertex[7]) +  
		tetrahedronVolume(vertex[1], vertex[3], vertex[4], vertex[6]);   
}


void Hexahedron::calculateCenter()
{
	double V1 =	tetrahedronVolume(vertex[0], vertex[1], vertex[3], vertex[4]);
	V3d c1 =	tetrahedronCenter(vertex[0], vertex[1], vertex[3], vertex[4]);

	double V2 =	tetrahedronVolume(vertex[1], vertex[2], vertex[3], vertex[6]);
	V3d c2 =	tetrahedronCenter(vertex[1], vertex[2], vertex[3], vertex[6]);

	double V3 =	tetrahedronVolume(vertex[1], vertex[4], vertex[5], vertex[6]);
	V3d c3 =	tetrahedronCenter(vertex[1], vertex[4], vertex[5], vertex[6]);

	double V4 =	tetrahedronVolume(vertex[3], vertex[4], vertex[6], vertex[7]);
	V3d c4 =	tetrahedronCenter(vertex[3], vertex[4], vertex[6], vertex[7]);

	double V5 =	tetrahedronVolume(vertex[1], vertex[3], vertex[4], vertex[6]);
	V3d c5 =	tetrahedronCenter(vertex[1], vertex[3], vertex[4], vertex[6]);

	center = (c1 * V1 + c2 * V2 + c3 * V3 + c4 * V4 + c5 * V5) / (V1 + V2 + V3 + V4 + V5);
}


Hexahedron::Hexahedron()
{
	numberOfVertex = 8;
	numberOfEdges = 12;
	numberOfFacets = 6;
}
