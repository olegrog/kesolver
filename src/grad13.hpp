#ifndef _GRAD13_HPP_
#define _GRAD13_HPP_

#include <stdexcept>

#include "ximesh.hpp"
#include "ximesh_mixture.hpp"
#include "ximesh_mix.hpp"
#include "ximesh_rot.hpp"
#include "distribution_function.hpp"

class Grad13 {
    public:
        typedef std::vector<std::string> Data;

        template <typename PropertyTree, typename XiMeshType> 
        Grad13(const PropertyTree& data,
                const XiMeshType& ximesh);
        
        size_t size() const { return f.size(); }
        double operator[](size_t i) const { return f[i]; }

        const DistributionFunction& func() const { return f; }

    private:
        DistributionFunction f;
};


template <typename F, Symmetry symmetry, template <Symmetry> class XiMeshType>
void setGrad13Simple(
        F& f,
        const typename XiMeshType<symmetry>::Vm v,
        const double temp,
        const double n,
        const typename XiMeshType<symmetry>::Vd t,
        const typename XiMeshType<symmetry>::Vm p,
        const typename XiMeshType<symmetry>::Vm q,
        const XiMeshType<symmetry>& mesh
)
{
    for (size_t i = 0; i < mesh.size(); ++i) {
        typename XiMeshType<symmetry>::Vm c = mesh.p(i) - v;
        f[i] = std::exp(-0.5 * sqr(c) / temp) * mesh.vol(i);
        double corr = 1 + (dot(c*c,n*(vd2vm(t)-temp))+2*dot(tr(c),p))/(2*n*sqr(temp)) +
            dot(q,c)*(sqr(c)/5/temp-1)/(n*sqr(temp));
        if (corr > 0) f[i] *= corr;
        if (f[i] <= 1e-100) std::cerr << "###Neg f(x): " << f[i] << '\n';
    }
}

template <typename F, Symmetry symmetry>
void setGrad13(
        F& f,
        const typename XiMesh<symmetry>::Vm v,
        const double temp,
        const double n,
        const typename XiMesh<symmetry>::Vd t,
        const typename XiMesh<symmetry>::Vm p,
        const typename XiMesh<symmetry>::Vm q,
        const XiMesh<symmetry>& mesh
)
{
    setGrad13Simple(f, v, temp, n, t, p, q, mesh);
}

template <typename F, Symmetry symmetry>
void setGrad13(
        F& f,
        const typename XiMesh<symmetry>::Vm v,
        const double temp,
        const double n,
        const typename XiMesh<symmetry>::Vd t,
        const typename XiMesh<symmetry>::Vm p,
        const typename XiMesh<symmetry>::Vm q,
        const XiMeshRect<symmetry>& mesh
)
{
    setGrad13Simple(f, v, temp, n, t, p, q, mesh);
}

template <typename F, Symmetry symmetry, template <Symmetry> class XiMesh>
void setGrad13Mixture(
        F& f,
        const typename XiMesh<symmetry>::Vm v,
        const double temp,
        const double n,
        const typename XiMesh<symmetry>::Vd t,
        const typename XiMesh<symmetry>::Vm p,
        const typename XiMesh<symmetry>::Vm q,
        const XiMesh<symmetry>& mesh
)
{
    throw std::invalid_argument("Grad13 not implemented.");
}

template <typename F, Symmetry symmetry>
void setGrad13(
        F& f,
        const typename XiMesh<symmetry>::Vm v,
        const double temp,
        const double n,
        const typename XiMesh<symmetry>::Vd t,
        const typename XiMesh<symmetry>::Vm p,
        const typename XiMesh<symmetry>::Vm q,
        const XiMeshMixture<symmetry>& mesh
)
{
    setGrad13Mixture(f, v, temp, n, t, p, q, mesh);
}

template <typename F, Symmetry symmetry>
void setGrad13(
        F& f,
        const typename XiMesh<symmetry>::Vm v,
        const double temp,
        const double n,
        const typename XiMesh<symmetry>::Vd t,
        const typename XiMesh<symmetry>::Vm p,
        const typename XiMesh<symmetry>::Vm q,
        const XiMeshMix<symmetry>& mesh
)
{
    setGrad13Mixture(f, v, temp, n, t, p, q, mesh);
}

template <typename F, Symmetry symmetry>
void setGrad13(
        F& f,
        const typename XiMesh<symmetry>::Vm v,
        const double temp,
        const double n,
        const typename XiMesh<symmetry>::Vd t,
        const typename XiMesh<symmetry>::Vm p,
        const typename XiMesh<symmetry>::Vm q,
        const XiMeshRot<symmetry>& mesh
)
{
    throw std::invalid_argument("Grad13 not implemented.");
}

template <typename PropertyTree, typename XiMeshType> 
Grad13::Grad13(
        const PropertyTree& data,
        const XiMeshType& ximesh
)
    : f(ximesh.size()) 
{
    typename XiMeshType::Vm v(strTo<typename XiMeshType::Vm>(data["u"].asString()));
    double temp = strTo<double>(data["T"].asString());
    double n = strTo<double>(data["n"].asString());
    typename XiMeshType::Vd t(strTo<typename XiMeshType::Vd>(data["t"].asString()));
    typename XiMeshType::Vm p(strTo<typename XiMeshType::Vm>(data["p"].asString()));
    typename XiMeshType::Vm q(strTo<typename XiMeshType::Vm>(data["q"].asString()));

    setGrad13(f, v, temp, n, t, p, q, ximesh);
    if (data.isMember("n")) {
        std::istringstream ss(data["n"].asString());
        equateN(f, ss, ximesh);
    }
}

#endif
