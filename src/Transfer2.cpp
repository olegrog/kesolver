#include <iostream>
#include "Transfer2.hpp"
#include "Constructors.hpp"

Transfer2::Transfer2(const MeshMpi& mesh)
{
	for (size_t i = 0; i < facets.size(); i++) 
		facets[i]->calculateDistance(spacemesh);
	for (size_t i = 0; i < spacemesh.size(); i++) 
		spacemesh[i]->inverseDD();
}

void Transfer2::move(const MeshMpi& mesh, const Gas& gas)
{
	data_exchanger.swap();

    for (size_t i = 0; i < spacemesh.size(); i++) {
        spacemesh[i]->prepareForNextStep();
    }
    for (size_t i = 0; i < facets.size(); i++) {
        facets[i]->findGradient(spacemesh, gas);
    }
    for (size_t i = 0; i < spacemesh.size(); i++) {
        spacemesh[i]->findGradient();
    }
    for (size_t i = 0; i < facets.size(); i++) {
        facets[i]->findPhi(spacemesh, gas);
    }
    for (size_t i = 0; i < facets.size(); i++) {
        facets[i]->transfer2(spacemesh, gas);
    }
    for (size_t i = 0; i < mypolys.size(); i++) 
        spacemesh[mypolys[i]]->f().equategf();

}

