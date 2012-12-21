#include <iostream>
#include <limits>

#include "Polygon.hpp"

void Polygon::calculateLength()
{
	Lmin = std::numeric_limits<double>::max();
	for (size_t i = 0; i < vertex.size(); ++i)
		for (size_t j = i+1; j < vertex.size(); ++j) {
			double l = norm(vertex[i] - vertex[j]);
			if (l < Lmin)
				Lmin = l;
		}
}

void Polygon::inverseDD() {
	dd = inverse(dd);
}

void Polygon::prepareForNextStep() {
    for (size_t i = 0; i < function.size(); ++i) {
        function.getFMin()[i]     = function.f()[i];
        function.getFMax()[i]     = function.f()[i];
        function.getGradient()[i] = 0.;
        function.getPhi()[i]      = 1.;
    }
}

void Polygon::findGradient() {
    for (size_t i = 0; i < function.size(); ++i)
        function.getGradient()[i] = mult(dd, function.getGradient()[i]);
}

