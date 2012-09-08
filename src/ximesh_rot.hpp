#ifndef _XIMESH_ROT_HPP_
#define _XIMESH_ROT_HPP_

#include <vector>
#include <iostream>

#include "symmetry.hpp"
#include "v.hpp"
#include "auxiliary.hpp"
#include "ximesh.hpp"

template <Symmetry symmetry>
struct RotVi {
    typedef typename SymmetryTrait<symmetry>::Vi Vi;
    RotVi() {}
    RotVi(const Vi vi, const int i) :
            vi(vi), i(i) {}
    operator Vi() const { return vi; }

    Vi vi;
    int i;
};

template <Symmetry symmetry>
inline RotVi<symmetry> operator+(const RotVi<symmetry> vi,
                                 const typename RotVi<symmetry>::Vi x)
{
    return RotVi<symmetry>(vi.vi + x, vi.i);
}

template <Symmetry symmetry>
inline RotVi<symmetry> operator-(const RotVi<symmetry> vi,
                                 const typename RotVi<symmetry>::Vi x)
{
    return RotVi<symmetry>(vi.vi - x, vi.i);
}

template <Symmetry symmetry>
struct RotVd {
    typedef typename SymmetryTrait<symmetry>::Vd Vd;
    RotVd() {}
    RotVd(const Vd v, const int i) :
            v(v), i(i) {}
    operator Vd() const { return v; }

    Vd v;
    int i;
};

template <Symmetry symmetry>
class XiMeshRot {
    public:
        XiMeshRot(int rad, double cut);

        size_t size() const { return xis.size(); }

        int radius() const { return rad;  }
        double a() const   { return a_;   }
        double cut() const { return cut_; }
        double vol() const { return vol_; }
        double vol(const int i) const;

        typedef typename SymmetryTrait<symmetry>::Vm Vm;
        typedef typename SymmetryTrait<symmetry>::Vd Vd;
        typedef RotVi<symmetry> Vi;
        typedef RotVd<symmetry> Vx;

        const Vd operator[](const int i) const 
                { return xis[i]; }

        const Vm p(const int i) const;

        double   e(const int i) const
                { return es[i]; }
        double erot(const int i) const
                { return erots[i]; }

        int operator()(const Vi i) const;

        const Vi i2vi(const int i) const
                { return vis[i]; }

        double i2q(const int i) const
                { return a_ * (i + 0.5); }

        const Vx i2xi(const Vi i) const
                { return Vx(::i2xi(i.vi, rad, a_), i.i); }

        const Vi xi2i(const Vx xi) const
                { return Vi(::xi2i(xi.v, rad, a_), xi.i); }

    private:
        int rad, vsize;
        double cut_, a_, vol_;
        std::vector<Vd> xis;
        std::vector<Vi> vis;
        std::vector<double> es, erots, qs;
        std::vector<int> xyzmap;

        int flatten(const Vi i) const
                { return i.i * vsize + ::flatten(i.vi, rad); }
};

template <>
inline const XiMeshRot<Cartesian>::Vm XiMeshRot<Cartesian>::p(const int i) const {
    return xis[i];
}

template <>
inline const XiMeshRot<Cylindrical>::Vm XiMeshRot<Cylindrical>::p(const int i) const {
    return xis[i][0];
}

template <>
inline double XiMeshRot<Cartesian>::vol(const int i) const {
    return qs[i];
}

template <>
inline double XiMeshRot<Cylindrical>::vol(const int i) const {
    return qs[i] * xis[i][1];
}

template <>
inline int XiMeshRot<Cartesian>::operator()(const Vi i) const {
    return ((i.vi >= 0) && (i.vi < V3i(2*rad)) && 
            (i.i  >= 0) && (i.i < rad)) ? xyzmap[flatten(i)] : -1;
}

template <>
inline int XiMeshRot<Cylindrical>::operator()(Vi i) const {
    if (i.vi[1] < 0) i.vi[1] = - i.vi[1] - 1;
    return ((i.vi >= 0) && (i.vi < V2i(2*rad, rad)) &&
            (i.i  >= 0) && (i.i < rad)) ? xyzmap[flatten(i)] : -1;
}

template <>
inline XiMeshRot<Cartesian>::XiMeshRot(int rad, double cut) : 
        rad(rad), vsize(cube(2*rad)), cut_(cut), a_(cut/rad) , vol_(cube(a_))
{
    std::cout << "vol_ = " << vol_ << std::endl;
    xyzmap.resize(rad * vsize);
    int i = 0;
    for (int k = 0; k < rad; ++k)
        for (int i1 = 0; i1 < 2*rad; ++i1)
            for (int i2 = 0; i2 < 2*rad; ++i2)
                for (int i3 = 0; i3 < 2*rad; ++i3) {
                    Vi xi(V3i(i1, i2, i3), k);
                    V3d v       = i2xi(xi);
                    double q    = i2q(k);
                    double erot = sqr(q);
                    double e    = sqr(v) + erot;
                    int j = flatten(xi);
                    if (e < sqr(cut_)) {
                        xis.push_back(v);
                        vis.push_back(xi);
                        es.push_back(e);
                        erots.push_back(erot);
                        qs.push_back(q);
                        xyzmap[j] = i;
                        ++i;
                    }
                    else {
                        xyzmap[j] = -1;
                    }
                }
}

template <>
inline XiMeshRot<Cylindrical>::XiMeshRot(int rad, double cut) : 
        rad(rad), vsize(2*sqr(rad)), cut_(cut), a_(cut/rad), vol_(2*M_PI*sqr(a_))
{
    xyzmap.resize(rad * vsize);
    int i = 0;
    for (int k = 0; k < rad; ++k)
        for (int i1 = 0; i1 < 2*rad; ++i1)
            for (int i2 = 0; i2 < rad; ++i2) {
                Vi xi(V2i(i1, i2), k);
                V2d v = i2xi(xi);
                double q    = i2q(k);
                double erot = sqr(q);
                double e    = sqr(v) + erot;
                int j = flatten(xi);
                if (e < sqr(cut_)) {
                    xis.push_back(v);
                    vis.push_back(xi);
                    es.push_back(e);
                    erots.push_back(erot);
                    qs.push_back(q);
                    xyzmap[j] = i;
                    ++i;
                }
                else {
                    xyzmap[j] = -1;
                }
            }
}

#endif
