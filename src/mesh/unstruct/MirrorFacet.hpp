#ifndef MIRRORFACET_H
#define MIRRORFACET_H

#include "PhysicalFacet.hpp"

class MirrorFacet : public PhysicalFacet {
	public:
		MirrorFacet()  {}

        void init(const PropertyTree& tree, const Gas& gas) {
			symm = static_cast<Axis>(tree["axis"].asInt());
        }

		void doTransfer(std::vector<Polygon*>& spacemesh, const Gas& gas);
		void doTransfer2(std::vector<Polygon*>& spacemesh, const Gas& gas);
		void calculateDistance(std::vector<Polygon*>& spacemesh);
		void doFindGradient(const std::vector<Polygon*>& spacemesh,
				const Gas& gas);

		int order() const { return 0; }
    private:
		Axis symm;

};

#endif /*MIRRORFACET_H*/
