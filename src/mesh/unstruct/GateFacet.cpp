#include "GateFacet.hpp"

#include "FacetFactory.hpp"

REGISTER_FACET(GateFacet, "gate")              

void GateFacet::doTransfer(std::vector<Polygon*>& spacemesh, const Gas& gas)
{
    size_t size       = spacemesh[polygon[0]]->f().size();

    SpeedFunction& f1 = spacemesh[polygon[0]]->f();
    DistributionFunction& f1_in     = f1.f();
    DistributionFunction& f2_in     = f1.g();

    SpeedFunction& f2 = spacemesh[polygon[1]]->f();
    DistributionFunction& f1_out    = f2.f();
    DistributionFunction& f2_out    = f2.g();

    const std::vector<V3d>& vels = gas.vel();

    for (size_t i = 0; i < size; ++i) {

        double sppr = dot(vels[i], n);
        if(sppr < 0.0) {
            f2_in[i]  += f1_in[i]*sppr*mult_in;
            f2_out[i] -= f1_in[i]*sppr*mult_out;
        }
        else{
            f2_in[i] += f1_out[i]*sppr*mult_in;
            f2_out[i] -= f1_out[i]*sppr*mult_out;
        }
    }
}

void GateFacet::doTransfer2(std::vector<Polygon*>& spacemesh, const Gas& gas)
{
    size_t size       = spacemesh[polygon[0]]->f().size();

    SpeedFunction& f1 = spacemesh[polygon[0]]->f();
    DistributionFunction&        f1_in     = f1.f();
    DistributionFunction&        f2_in     = f1.g();
    const DistributionFunction3& df_in     = f1.getGradient();
    const DistributionFunction&  phi_in    = f1.getPhi();

    SpeedFunction& f2 = spacemesh[polygon[1]]->f();
    DistributionFunction&        f1_out    = f2.f();
    DistributionFunction&        f2_out    = f2.g();
    const DistributionFunction3& df_out    = f2.getGradient();
    const DistributionFunction&  phi_out   = f2.getPhi();

    const std::vector<V3d>& vels = gas.vel();

    for (size_t i = 0; i < size; ++i) {
        double sppr = dot(vels[i], n);
        if (sppr < 0.0) {
            double dd = - 0.5 * dt * dot(vels[i], df_in[i]);
//            double dd = 0.0;
            double d = (f1_in[i] + (dot(df_in[i], d_in) + dd)*phi_in[i])*sppr;
            f2_in[i] += d * mult_in;
            f2_out[i] -= d * mult_out; 
            if ((f2_in[i] <= 0) || (f2_out[i] <= 0)) {
                std::cout << "GateFacet1!" << std::endl;
                std::cout << f1_in[i] << ' ' << f1_out[i] << std::endl;
                std::cout << f1_out[i] + (dot(df_out[i], d_out) + dd)*phi_out[i] << std::endl;
                std::cout << f2_in[i] << ' ' << f2_out[i] << std::endl;
                std::cout << phi_out[i] << std::endl;
                std::cout << f2.getFMax()[i] << ' ' << f2.getFMin()[i] << std::endl;
                std::cout << polygon[0] << ' ' << polygon[1] << std::endl;
                std::cout << i << std::endl;
                std::cout << vels[i] << std::endl;
                //exit(-1);
            }
        }
        else {
            double dd = - 0.5 * dt * dot(vels[i], df_out[i]);
//            double dd = 0.0;
            double d = (f1_out[i] + (dot(df_out[i], d_out) + dd)*phi_out[i])*sppr;
            f2_in[i] += d * mult_in; 
            f2_out[i] -= d * mult_out;
            if ((f2_in[i] <= 0) || (f2_out[i] <= 0)) {
                std::cout << "GateFacet2!" << std::endl;
                std::cout << f1_in[i] << ' ' << f1_out[i] << std::endl;
                std::cout << f1_out[i] + (dot(df_out[i], d_out) + dd)*phi_out[i] << std::endl;
                std::cout << f2_in[i] << ' ' << f2_out[i] << std::endl;
                std::cout << phi_out[i] << std::endl;
                std::cout << f2.getFMax()[i] << ' ' << f2.getFMin()[i] << std::endl;
                std::cout << polygon[0] << ' ' << polygon[1] << std::endl;
                std::cout << i << std::endl;
                std::cout << vels[i] << std::endl;
                //exit(-1);
            }
        }
    }
}

void GateFacet::calculateDistance(std::vector<Polygon*>& spacemesh) 
{
    Polygon* p_in  = spacemesh[polygon[0]];
    Polygon* p_out = spacemesh[polygon[1]];

    V3d d = p_out->getCenter() - p_in->getCenter();
    d_in  = getCenter() - p_in->getCenter();
    d_out = getCenter() - p_out->getCenter();

//  std::cout << "ds = " << d << std::endl;

    T3d dd = prod(d, d);
    p_in->getDD() += dd;
    p_out->getDD() += dd;
}

void GateFacet::doFindGradient(const std::vector<Polygon*>& spacemesh,
		const Gas& gas)
{
    size_t size       = spacemesh[polygon[0]]->f().size();

    SpeedFunction& f1 = spacemesh[polygon[0]]->f();
    const DistributionFunction& f_in = f1.f();
    DistributionFunction3& df_in     = f1.getGradient();
    DistributionFunction& fmin_in    = f1.getFMin();
    DistributionFunction& fmax_in    = f1.getFMax();

    SpeedFunction& f2 = spacemesh[polygon[1]]->f();
    const DistributionFunction& f_out = f2.f();
    DistributionFunction3& df_out     = f2.getGradient();
    DistributionFunction& fmin_out    = f2.getFMin();
    DistributionFunction& fmax_out    = f2.getFMax();

    V3d d_in_out = d_in - d_out;
    
    for (size_t i = 0; i < size; ++i) {
        V3d b = (f_out[i] - f_in[i]) * d_in_out;
        df_in[i] += b;
        df_out[i] += b;

        if (fmin_in[i] > f_out[i])
            fmin_in[i] = f_out[i];
        if (fmax_in[i] < f_out[i])
            fmax_in[i] = f_out[i];
        if (fmin_out[i] > f_in[i])
            fmin_out[i] = f_in[i];
        if (fmax_out[i] < f_in[i])
            fmax_out[i] = f_in[i];

    }
}

