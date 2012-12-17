#ifndef TRANSFER_H
#define TRANSFER_H

#include <vector>

#include "property_tree/property_tree.hpp"

#include "Gas.hpp"
#include "PhysicalFacet.hpp"
#include "DataExchanger.hpp"

class Transfer {
	protected:
		DataExchanger data_exchanger;
	public:
		virtual void init(const PropertyTree& tree, const Gas& gas, 
				const std::vector <PhysicalFacet*>& facets,	
				std::vector<Polygon*>& spacemesh, const std::vector<int>& mypolys, int rank) = 0;
		virtual	void move(const std::vector <PhysicalFacet*>& facets, 
				std::vector<Polygon*>& spacemesh, const std::vector<int>& mypolys, 
				const Gas& gas) = 0;
};

#endif /*TRANSFER_H*/
