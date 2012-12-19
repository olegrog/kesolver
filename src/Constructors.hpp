#ifndef _CONSTRUCTORS_HPP_
#define _CONSTRUCTORS_HPP_

#include <vector>

#include "property_tree/property_tree.hpp"

#include "Gas.hpp"

#include "mesh/unstruct/Polygon.hpp"

#include "Integral.hpp"

void GasConstructor(const PropertyTree& tree, Gas** gas_pp);

void GivePolygonMemoryAndInit(const PropertyTree& tree, const Gas& gas, Polygon* polygon); 

Integral IntegralConstructor(const PropertyTree& tree);

#endif // _CONSTRUCTORS_HPP_

