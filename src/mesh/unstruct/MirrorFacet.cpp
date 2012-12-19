#include "MirrorFacet.hpp"

#include "FacetFactory.hpp"

REGISTER_FACET(MirrorFacet, "mirror")              

void MirrorFacet::calculateDistance(std::vector<Polygon*>& spacemesh)
{
	Polygon* p_in = spacemesh[polygon[0]];

	d_in = 2.*(getCenter() - p_in->getCenter());

//	std::cout << "d_in (mirror) = " << d_in << std::endl;

	T3d dd = prod(d_in, d_in);
	p_in->getDD() += dd;
}

void MirrorFacet::doFindGradient(const std::vector<Polygon*>& spacemesh,
		const Gas& gas)
{
    size_t size       = spacemesh[polygon[0]]->f().size();

    SpeedFunction& f1 = spacemesh[polygon[0]]->f();
    const DistributionFunction& f_in      = f1.f();
    DistributionFunction3&      df_in     = f1.getGradient();
    DistributionFunction&       fmin_in   = f1.getFMin();
    DistributionFunction&       fmax_in   = f1.getFMax();

    for (size_t i = 0; i < size; ++i) {
		V3d b = (f_in[gas.mirror(i, symm)] - f_in[i]) * d_in;
		df_in[i] += b;

		if (fmin_in[i] > f_in[gas.mirror(i, symm)])
			fmin_in[i] = f_in[gas.mirror(i, symm)];
		if (fmax_in[i] < f_in[gas.mirror(i, symm)])
			fmax_in[i] = f_in[gas.mirror(i, symm)];
	}
}

void MirrorFacet::doTransfer(std::vector<Polygon*>& spacemesh, const Gas& gas)
{
    size_t size       = spacemesh[polygon[0]]->f().size();

    SpeedFunction& f1 = spacemesh[polygon[0]]->f();
    DistributionFunction& f1_in     = f1.f();
    DistributionFunction& f2_in     = f1.g();

	for (size_t i = 0; i < size; ++i) {
		double sppr = gas.dot(i, n);
		if (sppr < 0.0)
			f2_in[i] += f1_in[i]*sppr*mult_in;
		else
			f2_in[i] += f1_in[gas.mirror(i, symm)]*sppr*mult_in;

        if (f2_in[i] != f2_in[i]) {
            std::cout << "Nan Mirr\n";
            std::cout << f1_in[gas.mirror(i, symm)] << ' ' << f1_in[i] << std::endl;
            std::cout << sppr << std::endl;
            std::cout << mult_in << std::endl;
            exit(-1);
        }
	}
}

void MirrorFacet::doTransfer2(std::vector<Polygon*>& spacemesh, const Gas& gas)
{
    size_t size       = spacemesh[polygon[0]]->f().size();

    SpeedFunction& f1 = spacemesh[polygon[0]]->f();
    DistributionFunction&  f1_in       = f1.f();
    DistributionFunction&  f2_in       = f1.g();
    const DistributionFunction3& df_in = f1.getGradient();
    const DistributionFunction& phi_in = f1.getPhi();

	V3d d_in = getCenter() - spacemesh[polygon[0]]->getCenter();

    for (size_t i = 0; i < size; ++i) {
		double sppr = gas.dot(i, n);
		if (sppr < 0.0) {
            double dd = - 0.5 * dt * gas.dot(i, df_in[i]);
			f2_in[i] += (f1_in[i] + (dot(df_in[i], d_in)+dd)*phi_in[i])*sppr*mult_in;
        }
		else {
            double dd = - 0.5 * dt * gas.dot(i, df_in[gas.mirror(i, symm)]);
			f2_in[i] += (f1_in[gas.mirror(i, symm)] + 
					(dot(df_in[gas.mirror(i, symm)], d_in) + dd) * 
					    phi_in[gas.mirror(i, symm)]) * sppr * mult_in;
        }
	}
}

