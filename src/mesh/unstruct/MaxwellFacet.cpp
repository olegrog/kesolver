#include "MaxwellFacet.hpp"
#include "auxiliary.hpp"

#include "FacetFactory.hpp"

REGISTER_FACET(MaxwellFacet, "maxwell")              

void MaxwellFacet::calculateDistance(std::vector<Polygon*>& spacemesh)
{
	Polygon* p_in = spacemesh[polygon[0]];

	d_in = getCenter() - p_in->getCenter();

	T3d dd = prod(d_in, d_in);
	p_in->getDD() += dd;
}

void MaxwellFacet::doFindGradient(const std::vector<Polygon*>& spacemesh,
		const Gas& gas)
{
    size_t size       = spacemesh[polygon[0]]->f().size();

    SpeedFunction& f1 = spacemesh[polygon[0]]->f();
    const DistributionFunction& f_in = f1.f();
    DistributionFunction3& df_in     = f1.getGradient();
    DistributionFunction& fmin_in    = f1.getFMin();
    DistributionFunction& fmax_in    = f1.getFMax();

    for (size_t i = 0; i < size; ++i) {
		V3d b = (f[i] - f_in[i]) * d_in;
		df_in[i] += b;

		if (fmin_in[i] > f[i])
			fmin_in[i] = f[i];
		if (fmax_in[i] < f[i])
			fmax_in[i] = f[i];
	}
}

void MaxwellFacet::doTransfer(std::vector<Polygon*>& spacemesh, const Gas& gas)
{
    size_t size       = spacemesh[polygon[0]]->f().size();

    SpeedFunction& f1 = spacemesh[polygon[0]]->f();
    DistributionFunction& f1_in     = f1.f();
    DistributionFunction& f2_in     = f1.g();

	for (size_t i = 0; i < size; ++i) {
		double xin = gas.dot(i, n);
		if ( xin > 0. )
			f2_in[i] += f[i] * xin * mult_in;
		else
			f2_in[i] += f1_in[i] * xin * mult_in;
   }
}

void MaxwellFacet::doTransfer2(std::vector<Polygon*>& spacemesh, const Gas& gas)
{
    size_t size       = spacemesh[polygon[0]]->f().size();

    SpeedFunction& f1 = spacemesh[polygon[0]]->f();
    DistributionFunction& f1_in     = f1.f();
    DistributionFunction& f2_in     = f1.g();
    const DistributionFunction3& df = f1.getGradient();
    const DistributionFunction& phi = f1.getPhi();

	V3d dis = getCenter() - spacemesh[polygon[0]]->getCenter();

    for (size_t i = 0; i < size; ++i) {
		double xin = gas.dot(i, n);
		if ( xin > 0. )
			f2_in[i] += f[i] * xin * mult_in;
		else {
            double dd = - 0.5 * dt * gas.dot(i, df[i]);
			f2_in[i] += (f1_in[i] + (dot(df[i], dis)+dd)*phi[i]) * xin * mult_in;
        }
   }
}

