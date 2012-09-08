#ifndef _MAXWELL_HPP_
#define _MAXWELL_HPP_

#include "ximesh.hpp"
#include "ximesh_mixture.hpp"
#include "ximesh_mix.hpp"
#include "ximesh_rot.hpp"
#include "distribution_function.hpp"

class Maxwell {
    public:
        typedef std::vector<std::string> Data;

        template <typename XiMeshType> 
        Maxwell(const Data& data,
                const XiMeshType& ximesh);
        
        size_t size() const { return f.size(); }
        double operator[](size_t i) const { return f[i]; }

        const DistributionFunction& func() const { return f; }

    private:
        DistributionFunction f;
};


template <typename F, Symmetry symmetry>
void setMaxwell(F& f, const typename XiMesh<symmetry>::Vm v, const double temp,
                const XiMesh<symmetry>& mesh)
{
    const typename XiMesh<symmetry>::Vd u = vm2vd(v);
    for (size_t i = 0; i < mesh.size(); ++i)
        f[i] = std::exp( - 0.5 * sqr(mesh[i] - u) / temp) * mesh.vol(i);
}

template <typename F, Symmetry symmetry, template <Symmetry symmetry> class XiMeshType>
void setMaxwellMixture(F& f, const typename XiMeshType<symmetry>::Vm v, const double temp,
                const XiMeshType<symmetry>& mesh)
{
    typedef typename XiMesh<symmetry>::Vd Vd;
    const Vd u = vm2vd(v);
    for (size_t i = 0; i < mesh.size(); ++i) {
        const double m   = mesh.m(i);
        const Vd     xi  = mesh[i];
        const double vol = mesh.vol(i);
        f[i] = std::pow(m, -3./2.) * std::exp( - 0.5 * m * sqr(xi-u) / temp ) * vol;
    }
}

template <typename F, Symmetry symmetry>
void setMaxwell(F& f, const typename XiMesh<symmetry>::Vm v, const double temp,
                const XiMeshMixture<symmetry>& mesh)
{
    setMaxwellMixture(f, v, temp, mesh);
}

template <typename F, Symmetry symmetry>
void setMaxwell(F& f, const typename XiMesh<symmetry>::Vm v, const double temp,
                const XiMeshMix<symmetry>& mesh)
{
    setMaxwellMixture(f, v, temp, mesh);
}

template <typename F, Symmetry symmetry>
void setMaxwell(F& f, const typename XiMeshRot<symmetry>::Vm v, const double temp,
                const XiMeshRot<symmetry>& mesh)
{
    const typename XiMesh<symmetry>::Vd u = vm2vd(v);
    for (size_t i = 0; i < mesh.size(); ++i)
        f[i] = std::exp( - 0.5 * ( sqr(mesh[i] - u) + mesh.erot(i) ) / temp) * mesh.vol(i);
}

template <typename XiMeshType> 
Maxwell::Maxwell(const Data& data,
                 const XiMeshType& ximesh)
    : f(ximesh.size()) 
{
    std::string str;
    for (size_t i = 0; i < data.size(); ++i)
        str += data[i] + ' ';
	std::istringstream ss(str);
    typename XiMeshType::Vm v;
    double temp;
	ss >> v >> temp;

    setMaxwell(f, v, temp, ximesh);
    equateN(f, ss, ximesh);
}
    

#endif
