#include <iostream>
#include <vector>

#include "auxiliary.hpp"
#include "AlphaFacet.hpp"

#include "FacetFactory.hpp"

REGISTER_FACET(AlphaFacet, "alpha")              

void AlphaFacet::calculateDistance(std::vector<Polygon*>& spacemesh)
{
    Polygon* p_in = spacemesh[polygon[0]];

    d_in = getCenter() - p_in->getCenter();

    T3d dd = prod(d_in, d_in);

    idd = inverse(p_in->getDD());

    p_in->getDD() += dd;
}

void AlphaFacet::doFindGradient(const std::vector<Polygon*>& spacemesh,
        const Gas& gas)
{
    size_t size       = spacemesh[polygon[0]]->f().size();

    SpeedFunction& f1 = spacemesh[polygon[0]]->f();
    const DistributionFunction& f_in = f1.f();
    DistributionFunction3& df_in     = f1.getGradient();
    DistributionFunction& fmin_in    = f1.getFMin();
    DistributionFunction& fmax_in    = f1.getFMax();

    for (size_t i = 0; i < size; i++) {
        V3d df = mult(idd, df_in[i]);
        f_in2[i] = std::max(0., f_in[i] + dot(df, d_in));
    }

    gas.equateStreams(f, f_in2, n);

    for (size_t i = 0; i < size; ++i) {
        double sppr = gas.dot(i, n);
        double g1 = (sppr > 0.) ? f[i] : f_in2[i];
        double g2 = (f_in[gas.mirror(i, symm)] + f_in[i]) / 2;
        double g = alpha * g1 + (1 - alpha) * g2;
        V3d b = (g - f_in[i]) * d_in;

        df_in[i] += b;

        if (fmin_in[i] > g)
            fmin_in[i] = g;
        if (fmax_in[i] < g)
            fmax_in[i] = g;
    }
}

void AlphaFacet::doTransfer(std::vector<Polygon*>& spacemesh,
        const Gas& gas)
{
    size_t size       = spacemesh[polygon[0]]->f().size();

    SpeedFunction& f1 = spacemesh[polygon[0]]->f();
    DistributionFunction& f1_in     = f1.f();
    DistributionFunction& f2_in     = f1.g();

    gas.equateStreams(f, f1_in, n);

    for (size_t i = 0; i < size; ++i) {
        double sppr = gas.dot(i, n);

        if (sppr > 0.) {
            f2_in[i] += (alpha  * f[i] +
                    (1 - alpha) * f1_in[gas.mirror(i, symm)]) * sppr * mult_in;
        }
        else {
            f2_in[i] += f1_in[i] * sppr * mult_in;
        }
    }

}

void AlphaFacet::doTransfer2(std::vector<Polygon*>& spacemesh,
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
        if (sppr > 0.) {
            double dd = - 0.5 * dt * gas.dot(i, df_in[gas.mirror(i, symm)]);
            f2_in[i] += ( alpha  * f[i] +
                     (1 - alpha) * (f1_in[gas.mirror(i, symm)] + 
                    (dot(df_in[gas.mirror(i, symm)], d_in) + dd) * 
                        phi_in[gas.mirror(i, symm)]))
                     * sppr * mult_in;
        }
        else {
            double dd = - 0.5 * dt * gas.dot(i, df_in[i]);
            f2_in[i] += ( alpha  * f_in2[i] +
                     (1 - alpha) * (f1_in[i] + (dot(df_in[i], d_in)+dd)*phi_in[i]))
                     * sppr * mult_in;
        }

    }
}

