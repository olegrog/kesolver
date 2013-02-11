#include <iostream>
#include <limits>
#include <algorithm>

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

    std::copy ( function.f().begin(), function.f().end(), function.getFMin().begin() );
    std::copy ( function.f().begin(), function.f().end(), function.getFMax().begin() );

    std::fill( function.getGradient().begin(), function.getGradient().end(), V3d(0., 0., 0.) );
    std::fill( function.getPhi().begin(), function.getPhi().end(), 0. );

}

void Polygon::findGradientAndPhi(const Gas& gas) {

    DistributionFunction&    func = function.f();
    DistributionFunction&    g    = function.g();
    DistributionFunction&    fmax = function.getFMax();
    DistributionFunction&    fmin = function.getFMin();
    DistributionFunction3&  dfunc = function.getGradient();

    const std::vector<V3d>& vels = gas.vel();

    for (size_t i = 0; i < function.size(); ++i) {

        dfunc[i] = mult(dd, dfunc[i]);

        double f = func[i];
        V3d df = dfunc[i];
        g[i] = - 0.5 * dot(vels[i], df);
        fmax[i] -= f;
        fmin[i] -= f;

    }
}

