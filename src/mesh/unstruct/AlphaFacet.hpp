#ifndef ALPHAFACET_H
#define ALPHAFACET_H

#include "property_tree/property_tree.hpp"

#include "Gas.hpp"
#include "PhysicalFacet.hpp"
#include "distribution_function.hpp"
#include "v.hpp"

class AlphaFacet : public PhysicalFacet {
	public:
		AlphaFacet() {}

        void init(const PropertyTree& tree, const Gas& gas)
        {
            alpha = tree["alpha"].asDouble();
			symm = static_cast<Axis>(tree["axis"].asInt());
            DistributionFunction fm = gas.maxwell(tree);
            copy(fm, f);
            copy(fm, f_in2);
        }

		void doTransfer(std::vector<Polygon*>& spacemesh, const Gas& gas);
		void doTransfer2(std::vector<Polygon*>& spacemesh, const Gas& gas);
		void calculateDistance(std::vector<Polygon*>& spacemesh);
		void doFindGradient(const std::vector<Polygon*>& spacemesh,
				const Gas& gas);
		int order() const { return 1; }

	private:
        double alpha;
		Axis symm;
		T3d idd;
        DistributionFunction f, f_in2;
};


#endif /*ALPHAFACET_H*/

