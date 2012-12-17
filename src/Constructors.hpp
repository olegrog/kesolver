#ifndef _CONSTRUCTORS_HPP_
#define _CONSTRUCTORS_HPP_

#include <vector>

#include "property_tree/property_tree.hpp"

#include "Gas.hpp"

#include "Polygon.hpp"

#include "PhysicalFacet.hpp"

#include "Integral.hpp"

void GasConstructor(const PropertyTree& tree, Gas** gas_pp);

void ElementsConstructor(const PropertyTree& tree,
                         std::vector<Polygon*>& polygons,
                         std::vector<PhysicalFacet*>& facets, 
                         const Gas& gas);

void MypolysConstructor(int rank, 
                        const std::vector<Polygon*>& polygons, 
                        std::vector<int>& mypolys);

double findTimeStep(const std::vector<Polygon*>& spacemesh, 
                    const Gas& gas, 
                    double curnt);

void GivePolygonMemoryAndInit(const PropertyTree& tree, const Gas& gas, Polygon* polygon); 

Integral IntegralConstructor(const PropertyTree& tree);

#endif // _CONSTRUCTORS_HPP_

