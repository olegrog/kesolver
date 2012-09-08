#ifndef WALLMAXWELLFACET_H
#define WALLMAXWELLFACET_H

#include "Gas.hpp"
#include "PhysicalFacet.hpp"
#include "distribution_function.hpp"
#include "v.hpp"

class WallMaxwellFacet : public PhysicalFacet {
	public:
        template <typename F>
		WallMaxwellFacet(const F& f_)
                { copy(f_, f); copy(f_, f_in2); }

		void doTransfer(std::vector<Polygon*>& spacemesh, const Gas& gas);
		void doTransfer2(std::vector<Polygon*>& spacemesh, const Gas& gas);
		void calculateDistance(std::vector<Polygon*>& spacemesh);
		void doFindGradient(const std::vector<Polygon*>& spacemesh,
				const Gas& gas);
		int order() const { return 1; }

	private:
		T3d idd;

        DistributionFunction f, f_in2;
};


#endif /*WALLMAXWELLFACET_H*/
