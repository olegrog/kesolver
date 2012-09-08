#ifndef _MAXWELLFACET_HPP_
#define _MAXWELLFACET_HPP_

#include "PhysicalFacet.hpp"
#include "distribution_function.hpp"
#include "v.hpp"

class MaxwellFacet : public PhysicalFacet {
	public:
        template <typename F>
		MaxwellFacet(const F& f_)
                { copy(f_, f); }

		void doTransfer(std::vector<Polygon*>& spacemesh, const Gas& gas);
		void doTransfer2(std::vector<Polygon*>& spacemesh, const Gas& gas);
		void calculateDistance(std::vector<Polygon*>& spacemesh);
		void doFindGradient(const std::vector<Polygon*>& spacemesh,
				const Gas& gas);
		int order() const { return 0; }

	private:
        DistributionFunction f;
        
};


#endif // _MAXWELLFACET_HPP_
