#ifndef TRANSFER2_HPP_
#define TRANSFER2_HPP_

#include "Transfer.hpp"

class Transfer2 : public Transfer {
	public:
		Transfer2(const std::vector <PhysicalFacet*>& facets, 
				std::vector<Polygon*>& spacemesh, const std::vector<int>& mypolys);
		void move(const std::vector <PhysicalFacet*>& facets,
				std::vector<Polygon*>& spacemesh, const std::vector<int>& mypolys,
				const Gas& gas);
		void init(const Loader& ldr, const Gas& gas, 
				const std::vector <PhysicalFacet*>& facets,	
				std::vector<Polygon*>& spacemesh, const std::vector<int>& mypolys, int rank);
};

#endif /*TRANSFER2_HPP_*/
