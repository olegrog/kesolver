#ifndef GATEFACET_H
#define GATEFACET_H

#include "PhysicalFacet.hpp"

class GateFacet : public PhysicalFacet {
	public:
		void calculateDistance(std::vector<Polygon*>& spacemesh);

		void doTransfer(std::vector<Polygon*>& spacemesh, const Gas& gas);
		void doTransfer2(std::vector<Polygon*>& spacemesh, const Gas& gas);

		void doFindGradient(const std::vector<Polygon*>& spacemesh,
				const Gas& gas);

		int order() const { return 0; }
};

#endif /*GATEFACET_H*/
