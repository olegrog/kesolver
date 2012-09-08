#ifndef _XIMESH_MIXTURE_HPP_
#define _XIMESH_MIXTURE_HPP_

#include <vector>
#include <algorithm>
#include <cassert>

#include "ximesh.hpp"

template <Symmetry symmetry>
struct MixtureVi {
    typedef typename SymmetryTrait<symmetry>::Vi Vi;
    MixtureVi() {}
    MixtureVi(const Vi vi, const int i) :
            vi(vi), i(i) {}
    operator Vi() const { return vi; }

    Vi vi;
    int i;
};

template <Symmetry symmetry>
inline MixtureVi<symmetry> operator+(const MixtureVi<symmetry> vi,
                                     const typename MixtureVi<symmetry>::Vi x)
{
    return MixtureVi<symmetry>(vi.vi + x, vi.i);
}

template <Symmetry symmetry>
inline MixtureVi<symmetry> operator-(const MixtureVi<symmetry> vi,
                                     const typename MixtureVi<symmetry>::Vi x)
{
    return MixtureVi<symmetry>(vi.vi - x, vi.i);
}

template <Symmetry symmetry>
struct MixtureVd {
    typedef typename SymmetryTrait<symmetry>::Vd Vd;
    MixtureVd() {}
    MixtureVd(const Vd v, const int i) :
            v(v), i(i) {}
    operator Vd() const { return v; }

    Vd v;
    int i;
};

template <Symmetry symmetry>
class XiMeshMixture {
    public:
        template <typename VecI, typename VecD>
        XiMeshMixture(double a, const VecI& rads, const VecD& ms);

        size_t size() const { return size_; }

        double vol()            const { return vol_; }
        double vol(const int i) const { return vols[i]; }

        typedef typename SymmetryTrait<symmetry>::Vm Vm;
        typedef typename SymmetryTrait<symmetry>::Vd Vd;
        typedef typename SymmetryTrait<symmetry>::Vj Vj;
        typedef MixtureVi<symmetry> Vi;
        typedef MixtureVd<symmetry> Vx;

        const Vd operator[](const int i) const 
                { return xis[i]; }

        const Vm p(const int i) const { return ps[i]; }
        double   e(const int i) const { return es[i]; }
        double   m(const int i) const { return ms[i]; }

        int operator()(const Vi i) const {
            const int j = ximeshes[i.i](i.vi);
            return (j >= 0) ? flatten(i.i, j) : -1;
        }

        const Vi i2vi(const int i) const { return vis[i]; }
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

        const std::vector<int>&    rads() const { return rads_; }
        const std::vector<double>& cuts() const { return cuts_; }
        const std::vector<double>& masses() const { return masses_; }

        int mirror(int i, const Axis a) const;

    private:
        XiMeshes ximeshes;
        std::vector<double> cuts_, masses_;
        std::vector<int> offsets, rads_;
        size_t size_;
        double vol_, cut_;
        std::vector<Vd> xis;
        std::vector<Vi> vis;
        std::vector<Vm> ps;
        std::vector<double> es, vols, ms;
        std::vector<int> cis;
        std::vector<Vj> mirr;

        int flatten(const int i_mix, const int i_xi) const {
            return offsets[i_mix] + i_xi;
        }

};

template <>
inline int XiMeshMixture<Cartesian>::mirror(const int i, const Axis a) const {
    return mirr[i][a];
}

template <>
inline int  XiMeshMixture<Cylindrical>::mirror(const int i, const Axis a) const {
    return mirr[i];
}


template <Symmetry symmetry> template <typename VecI, typename VecD>
XiMeshMixture<symmetry>::XiMeshMixture(double a, const VecI& rads, const VecD& masses) :
        masses_(masses), rads_(rads)
{
    assert(rads.size() == masses.size());

    ximeshes.reserve(rads.size());
    offsets.reserve(rads.size()+1);

    int offset = 0;
    for (size_t i = 0; i < rads.size(); ++i) 
    {
        const int rad = rads[i];
        const double mass = masses[i];
        std::cout << "rad, cut = " << rad << ' ' << a*rad << std::endl;
        XiMeshType ximesh(rad, a*rad);
        offsets.push_back(offset);
        offset += ximesh.size();
        ximeshes.push_back(ximesh);
        cuts_.push_back(ximesh.cut() / mass);
    }
    cut_  = *std::max_element(cuts_.begin(), cuts_.end());
    vol_  = ximeshes[0].vol();
    offsets.push_back(offset);
    size_ = offset;

    std::cout << "size = " << size_ << std::endl;

    xis.reserve(size_);
    vis.reserve(size_);
    ps.reserve(size_);
    es.reserve(size_);
    vols.reserve(size_);
    ms.reserve(size_);
    cis.reserve(size_);

    for (size_t j = 0; j < ximeshes.size(); ++j) {
        const XiMeshType& mesh = ximeshes[j];
        for (size_t i = 0; i < mesh.size(); ++i) {
             xis.push_back(mesh[i] / masses[j]);
             vis.push_back(Vi(mesh.i2vi(i), j));
              ps.push_back(mesh.p(i));
              es.push_back(mesh.e(i) / masses[j]);
            vols.push_back(mesh.vol(i));
              ms.push_back(masses[j]);
             cis.push_back(j);
            mirr.push_back(offsets[j] + mesh.mirror(i));
        }
    }
}

#endif
