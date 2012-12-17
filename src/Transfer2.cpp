#include <iostream>
#include "Transfer2.hpp"
#include "Constructors.hpp"

Transfer2::Transfer2(const std::vector <PhysicalFacet*>& facets, std::vector<Polygon*>& spacemesh, const std::vector<int>& mypolys)
{
	for (size_t i = 0; i < facets.size(); i++) 
		facets[i]->calculateDistance(spacemesh);
	for (size_t i = 0; i < spacemesh.size(); i++) 
		spacemesh[i]->inverseDD();
}

void Transfer2::move(const std::vector <PhysicalFacet*>& facets, std::vector<Polygon*>& spacemesh, const std::vector<int>& mypolys, const Gas& gas)
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

void Transfer2::init(const PropertyTree& tree, const Gas& gas, 
				const std::vector <PhysicalFacet*>& facets,	
				std::vector<Polygon*>& spacemesh, const std::vector<int>& mypolys, int rank)
{
	data_exchanger.init2(spacemesh, mypolys, rank);
	for (size_t i = 0; i < mypolys.size(); ++i) {
		GivePolygonMemoryAndInit(tree, gas, spacemesh[mypolys[i]]);
		spacemesh[mypolys[i]]->f().giveMemoryToGradient();
	}
	const std::vector<int>& toAllocPolygons = data_exchanger.getToAllocPolygons();
	for (size_t i = 0; i < toAllocPolygons.size(); ++i) {
		GivePolygonMemoryAndInit(tree, gas, spacemesh[toAllocPolygons[i]]);
		spacemesh[toAllocPolygons[i]]->f().giveMemoryToGradient();
	}

	for (size_t i = 0; i < facets.size(); ++i) {
		facets[i]->activation(spacemesh);
	}

	data_exchanger.mpiInit(spacemesh);
}

