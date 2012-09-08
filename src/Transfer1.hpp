#ifndef TRANSFER1_HPP_
#define TRANSFER1_HPP_

#include "Transfer.hpp"

class Transfer1 : public Transfer {
	public:
		void move(const std::vector <PhysicalFacet*>& facets,
				std::vector<Polygon*>& spacemesh, const std::vector<int>& mypolys,
				const Gas& gas);
		void init(const Loader& ldr, const Gas& gas, 
				const std::vector <PhysicalFacet*>& facets,	
				std::vector<Polygon*>& spacemesh, const std::vector<int>& mypolys, int rank);
};

#endif /*TRANSFER1_HPP_*/
