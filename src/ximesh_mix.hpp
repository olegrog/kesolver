#ifndef _XIMESH_MIX_HPP_
#define _XIMESH_MIX_HPP_

#include <vector>
#include <algorithm>

#include "ximesh.hpp"
#include "ximesh_mixture.hpp"

template <Symmetry symmetry>
class XiMeshMix {
    public:
        typedef typename SymmetryTrait<symmetry>::Vm Vm;
        typedef typename SymmetryTrait<symmetry>::Vd Vd;
        typedef MixtureVi<symmetry> Vi;
        typedef MixtureVd<symmetry> Vx;
        typedef typename SymmetryTrait<symmetry>::Vj Vj;

        template <typename VecD>
        XiMeshMix(const int rad, const double s2e, const VecD& ms,
                  const Vm v_, const VecD& adds);

        size_t size() const { return size_; }

        double vol()            const { return vol_; }
        double vol(const int i) const { return vols[i]; }
        double volLog(const int i) const { return vol_logs[i]; }

        double a(const int i) const { return as[i]; }

        const Vd operator[](const int i) const 
                { return xis[i]; }

        const std::vector<V3d>& vel() const 
                { return v3; }

        const Vm p(const int i) const { return ps[i]; }
        double   e(const int i) const { return es[i]; }
        double   m(const int i) const { return ms[i]; }

        int operator()(const Vi i) const {
            const int j = ximeshes[i.i](i.vi);
            return (j >= 0) ? flatten(i.i, j) : -1;
        }

        const Vi  i2vi(const int i) const { return vis[i]; }
        const int i2ci(const int i) const { return cis[i]; }

        const Vx i2xi(const Vi i) const {
            return Vx(ximeshes[i.i].i2xi(i.vi), i.i);
        }
        const Vi xi2i(const Vx xi) const {
            return Vi(ximeshes[xi.i].xi2i(xi.v), xi.i);
        }

        typedef XiMesh<symmetry> XiMeshType;
        typedef std::vector< XiMeshType > XiMeshes;

        const XiMeshes& xiMeshes() const { return ximeshes; }

        const int offset(const int i) const { return offsets[i]; }

        double cut() const { return cut_; }

        int   rad() const { return rad_; }
        const std::vector<double>& cuts() const { return cuts_; }
        const std::vector<double>& masses() const { return masses_; }

        int mirror(int i, const Axis a) const;

    private:
        XiMeshes ximeshes;
        std::vector<double> masses_, cuts_;
        std::vector<int> offsets;
        int size_, rad_;
        double vol_, cut_;
        std::vector<Vd> xis;
        std::vector<Vi> vis;
        std::vector<Vm> ps;
        std::vector<double> es, vols, ms, vol_logs, as;
        std::vector<int> cis;
        std::vector<Vj> mirr;
        std::vector<V3d> v3;

        int flatten(const int i_mix, const int i_xi) const {
            return offsets[i_mix] + i_xi;
        }

};

template <>
inline int XiMeshMix<Cartesian>::mirror(const int i, const Axis a) const {
    return mirr[i][a];
}

template <>
inline int  XiMeshMix<Cylindrical>::mirror(const int i, const Axis a) const {
    return mirr[i];
}

template <Symmetry symmetry> template <typename VecD>
XiMeshMix<symmetry>::XiMeshMix(const int rad, const double s2e, const VecD& masses, 
                               const Vm v_, const VecD& adds) :
        masses_(masses), rad_(rad)
{
    assert(adds.size() == masses.size());

    ximeshes.reserve(masses.size());
    offsets.reserve(masses.size()+1);
    cuts_.reserve(masses.size());

    int offset = 0;
    for (size_t i = 0; i < masses.size(); ++i) 
    {
        const double mass = masses[i];
        const double cut = std::sqrt(mass) * s2e + mass * adds[i];
        std::cout << "cut = " << cut << std::endl;
        const Vm v = v_ * mass;
        XiMeshType ximesh(rad, cut, v);
        offsets.push_back(offset);
        offset += ximesh.size();
        ximeshes.push_back(ximesh);
        cuts_.push_back(cut / mass);
    }
    cut_  = *std::max_element(cuts_.begin(), cuts_.end());
    vol_  = 1;
    offsets.push_back(offset);
    size_ = offset;

    std::cout << "size = " << size_ << std::endl;

    xis.reserve(size_);
    vis.reserve(size_);
    ps.reserve(size_);
    es.reserve(size_);
    vols.reserve(size_);
    vol_logs.reserve(size_);
    ms.reserve(size_);
    cis.reserve(size_);
    as.reserve(size_);
    v3.reserve(size_);

    for (size_t j = 0; j < ximeshes.size(); ++j) {
        const XiMeshType& mesh = ximeshes[j];
        for (size_t i = 0; i < mesh.size(); ++i) {
             xis.push_back(mesh[i] / masses[j]);
             vis.push_back(Vi(mesh.i2vi(i), j));
              ps.push_back(mesh.p(i));
              es.push_back(mesh.e(i) / masses[j]);
            double vol = mesh.vol(i) * mesh.vol();
            vols.push_back(vol);
        vol_logs.push_back(std::log(vol));
              ms.push_back(masses[j]);
             cis.push_back(j);
              as.push_back(mesh.a());
            mirr.push_back(offsets[j] + mesh.mirror(i));
              v3.push_back(mesh.vel()[i] / masses[j]);
        }
    }
}


#endif
