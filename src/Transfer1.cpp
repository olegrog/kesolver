#include <iostream>
#include "Transfer1.hpp"
#include "Constructors.hpp"

void Transfer1::move(const std::vector <PhysicalFacet*>& facets, std::vector<Polygon*>& spacemesh, const std::vector<int>& mypolys, const Gas& gas)
{
	data_exchanger.swap();

    for(size_t i = 0; i < facets.size(); ++i)
        facets[i]->transfer(spacemesh, gas);
    for(size_t i = 0; i < mypolys.size(); ++i) 
        spacemesh[mypolys[i]]->f().equategf();

}

void Transfer1::init(const Loader& ldr, const Gas& gas,
				const std::vector <PhysicalFacet*>& facets,	
				std::vector<Polygon*>& spacemesh, const std::vector<int>& mypolys, int rank)
{
	data_exchanger.init(spacemesh, mypolys, rank);
	for (size_t i = 0; i < mypolys.size(); ++i) {
		GivePolygonMemoryAndInit(ldr, gas, spacemesh[mypolys[i]]);
	}
	const std::vector<int>& toAllocPolygons = data_exchanger.getToAllocPolygons();
	for (size_t i = 0; i < toAllocPolygons.size(); ++i) {
		GivePolygonMemoryAndInit(ldr, gas, spacemesh[toAllocPolygons[i]]);
	}
	for (size_t i = 0; i < facets.size(); ++i) {
		facets[i]->activation(spacemesh);
	}
	data_exchanger.mpiInit(spacemesh);
}

