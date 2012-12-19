#ifndef INTEGRAL_H
#define INTEGRAL_H

#include <vector>
#include "mesh/unstruct/Polygon.hpp"
#include "Gas.hpp"
#include "section.hpp"

class Integral{

	public:
		void collide(double t, 
                     std::vector<Polygon*>& spacemesh, 
                     const std::vector<int>& mypolys, 
                     Gas& gas);

		Integral(int p_, int order_, 
                bool is_free_molecular,
                const SimpleSection* section);

	private:
		int p;
        int order;
        bool is_free_molecular;
        const SimpleSection* section;
};

#endif /*INTEGRAL_H*/
