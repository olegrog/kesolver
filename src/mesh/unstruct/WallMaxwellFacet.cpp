#include <iostream>
#include <vector>

#include "auxiliary.hpp"
#include "WallMaxwellFacet.hpp"

#include "FacetFactory.hpp"

REGISTER_FACET(WallMaxwellFacet, "diffusion")              

void WallMaxwellFacet::calculateDistance(std::vector<Polygon*>& spacemesh)
{
	Polygon* p_in = spacemesh[polygon[0]];

	d_in = getCenter() - p_in->getCenter();

//	std::cout << "d_in = " << d_in << std::endl;

	T3d dd = prod(d_in, d_in);

	idd = inverse(p_in->getDD());
//	std::cout << "idd = " << idd << std::endl;

	p_in->getDD() += dd;
}

void WallMaxwellFacet::doFindGradient(const std::vector<Polygon*>& spacemesh,
		const Gas& gas)
{
    size_t size       = spacemesh[polygon[0]]->f().size();

    SpeedFunction& f1 = spacemesh[polygon[0]]->f();
    const DistributionFunction& f_in = f1.f();
    DistributionFunction3& df_in     = f1.getGradient();
    DistributionFunction& fmin_in    = f1.getFMin();
    DistributionFunction& fmax_in    = f1.getFMax();

    gas.equateStreams(f, f_in, n);

    for (size_t i = 0; i < size; ++i) {
		double sppr = gas.dot(i, n);
		double g;
		if (sppr > 0.) 
			g = f[i];
		else {
			V3d df = mult(idd, df_in[i]);
			g = std::max(0., f_in[i] + dot(df, d_in));
		}

		V3d b = (g - f_in[i]) * d_in;
		df_in[i] += b;

		if (fmin_in[i] > g)
			fmin_in[i] = g;
		if (fmax_in[i] < g)
			fmax_in[i] = g;
	}

}

void WallMaxwellFacet::doTransfer(std::vector<Polygon*>& spacemesh,
		const Gas& gas)
{
    size_t size       = spacemesh[polygon[0]]->f().size();

    SpeedFunction& f1 = spacemesh[polygon[0]]->f();
    DistributionFunction& f1_in     = f1.f();
    DistributionFunction& f2_in     = f1.g();

    gas.equateStreams(f, f1_in, n);

	for (size_t i = 0; i < size; ++i) {
		double sppr = gas.dot(i, n);

		if (sppr > 0.) 
			f2_in[i] += f[i] * sppr * mult_in;
        else
			f2_in[i] += f1_in[i] * sppr * mult_in;
	}

}

void WallMaxwellFacet::doTransfer2(std::vector<Polygon*>& spacemesh,
		const Gas& gas)
{
    size_t size       = spacemesh[polygon[0]]->f().size();

    SpeedFunction& f1 = spacemesh[polygon[0]]->f();
    DistributionFunction& f1_in        = f1.f();
    DistributionFunction& f2_in        = f1.g();
    const DistributionFunction3& df_in = f1.getGradient();
    const DistributionFunction& phi_in = f1.getPhi();

	for (size_t i = 0; i < size; i++) {
        double dd = - 0.5 * dt * gas.dot(i, df_in[i]);
        f_in2[i] = f1_in[i] + (dot(df_in[i], d_in) + dd) * phi_in[i];
    }

    gas.equateStreams(f, f_in2, n);

    for (size_t i = 0; i < size; ++i) {
		double sppr = gas.dot(i, n);
		if (sppr > 0.) 
			f2_in[i] += f[i] * sppr * mult_in;
        else {
			double d = f_in2[i] * sppr;
			f2_in[i] += d * mult_in;
        }
	}

}

