#ifndef ALPHAFACET_H
#define ALPHAFACET_H

#include "Gas.hpp"
#include "PhysicalFacet.hpp"
#include "distribution_function.hpp"
#include "v.hpp"

class AlphaFacet : public PhysicalFacet {
	public:
        template <typename F>
		AlphaFacet(double a, Axis s, const F& f_) : 
            alpha(a), symm(s)
        {
            copy(f_, f); copy(f_, f_in2);
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
