#ifndef _CONSTRUCTORS_HPP_
#define _CONSTRUCTORS_HPP_

#include <vector>
#include "Loader.hpp"
#include "Polygon.hpp"
#include "PhysicalFacet.hpp"
#include "Gas.hpp"
#include "Integral.hpp"

void PolygonsConstructor(const Loader& loader, std::vector<Polygon*>& polygons);

void FacetsConstructor(const Loader& loader, const std::vector<Polygon*>& polygons,
		std::vector<PhysicalFacet*>& facets, 
        const Gas& gas,
        double time_step, 
		int rank);

void MypolysConstructor(int rank, const std::vector<Polygon*>& polygons, std::vector<int>& mypolys);

void GivePolygonMemoryAndInit(const Loader& loader, const Gas& gas, Polygon* polygon); 

double findTimeStep(const std::vector<Polygon*>& spacemesh, const Gas& gas, double curnt);

void GasConstructor(const Loader& loader, Gas** gas_pp);

Integral IntegralConstructor(const Loader& loader);

#endif // _CONSTRUCTORS_HPP_

