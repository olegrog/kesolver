#ifndef _XIMESH_HPP_
#define _XIMESH_HPP_

#include <vector>
#include <iostream>

#include "symmetry.hpp"
#include "v.hpp"
#include "auxiliary.hpp"
#include "axis.hpp"

inline const V3d i2xi(const V3i i, const int rad, const double a) {
    return (V3d(i)+0.5-rad) * a;
}

inline const V2d i2xi(const V2i i, const int rad, const double a) {
    return (V2d(i)+0.5-V2d(rad, 0)) * a;
}

inline const V3i xi2i(const V3d xi, const int rad, const double a) {
    return V3i(xi/a+rad+1000)-1000;
}

inline const V2i xi2i(const V2d xi, const int rad, const double a) {
    return V2i(xi/a+V2d(rad, 0)+1000)-1000;
}

inline int flatten(const V3i i, const int rad) {
    return dot(i, V3i(1, 2*rad, sqr(2*rad)));
}

inline int flatten(const V2i i, const int rad) {
    return dot(i, V2i(1, 2*rad));
}

template <Symmetry symmetry>
class XiMesh {
    public:
        typedef typename SymmetryTrait<symmetry>::Vm Vm;
        typedef typename SymmetryTrait<symmetry>::Vd Vd;
        typedef typename SymmetryTrait<symmetry>::Vi Vi;
        typedef typename SymmetryTrait<symmetry>::Vj Vj;
        typedef Vd Vx;

        XiMesh(int rad, double cut, const Vm v = 0.);

        XiMesh(const XiMesh<symmetry>& ximesh) :
            rad(ximesh.rad), cut_(ximesh.cut_),
            a_(ximesh.a_), vol_(ximesh.vol_), 
            xis(ximesh.xis), vis(ximesh.vis),
            xyzmap(ximesh.xyzmap), mirr(ximesh.mirr),
            v3(ximesh.v3)
        {
            std::cout << "XiMesh::CopyConstructor\n";
        }

        size_t size() const { return xis.size(); }

        int radius() const { return rad;  }
        double a() const   { return a_;   }
        double cut() const { return cut_ + max(shift); }
        double vol() const { return vol_; }
        double vol(const int i) const;

        const Vd operator[](const int i) const 
                { return xis[i]; }

        const std::vector<V3d>& vel() const {
            return v3;
        }

        const Vm p(const int i) const;
        double   e(const int i) const 
                { return sqr(operator[](i)); }

        int operator()(const Vi i) const;
        const Vi i2vi(const int i) const { return vis[i]; }

        const Vx i2xi(const Vi i) const
                { return ::i2xi(i, rad, a_) + shift; }
        const Vi xi2i(const Vx xi) const
                { return ::xi2i(xi - shift, rad, a_); }

        int mirror(int i, const Axis a) const;
        const Vj mirror(int i) const
                { return mirr[i]; }

    private:
        int rad;
        double cut_, a_, vol_;
        Vd shift;
        std::vector<Vd> xis;
        std::vector<Vi> vis;
        std::vector<int> xyzmap;
        std::vector<Vj> mirr;
        std::vector<V3d> v3;

        int flatten(const Vi i) const 
                { return ::flatten(i, rad); }
};

template <>
inline const XiMesh<Cartesian>::Vm XiMesh<Cartesian>::p(const int i) const {
    return xis[i];
}

template <>
inline const XiMesh<Cylindrical>::Vm XiMesh<Cylindrical>::p(const int i) const {
    return xis[i][0];
}

template <>
inline double XiMesh<Cartesian>::vol(const int i) const {
    return 1;
}

template <>
inline double XiMesh<Cylindrical>::vol(const int i) const {
    return xis[i][1];
}

template <>
inline int XiMesh<Cartesian>::operator()(const Vi i) const {
    return ((i >= 0) && (i < Vi(2*rad))) ? xyzmap[flatten(i)] : -1;
}

template <>
inline int XiMesh<Cylindrical>::operator()(Vi i) const {
    if (i[1] < 0) i[1] = - i[1] - 1;
    return ((i >= 0) && (i < Vi(2*rad, rad))) ? xyzmap[flatten(i)] : -1;
}

template <>
inline int XiMesh<Cartesian>::mirror(const int i, const Axis a) const {
    return mirr[i][a];
}

template <>
inline int  XiMesh<Cylindrical>::mirror(const int i, const Axis a) const {
    return mirr[i];
}


template <>
inline XiMesh<Cartesian>::XiMesh(int rad, double cut, const Vm v_ ) : 
        rad(rad), cut_(cut), a_(cut/rad) , vol_(cube(a_)), shift(vm2vd(v_))
{
    std::cout << "vol_ = " << vol_ << std::endl;
    xyzmap.resize(cube(2*rad));
    int i = 0;
    for (int i1 = 0; i1 < 2*rad; ++i1)
        for (int i2 = 0; i2 < 2*rad; ++i2)
            for (int i3 = 0; i3 < 2*rad; ++i3) {
                V3i xi(i1, i2, i3);
                V3d v = i2xi(xi) - shift;
                int j = flatten(xi);
                if (sqr(v) < sqr(cut_)) {
                    xis.push_back(v + shift);
                    vis.push_back(xi);
                    xyzmap[j] = i;
                    ++i;
                }
                else {
                    xyzmap[j] = -1;
                }
            }

    mirr.resize(vis.size());
    for (size_t i = 0; i < vis.size(); ++i) {
        V3i xi(vis[i]);
        mirr[i] = V3i(operator()(V3i(2*rad-1-xi[0], xi[1], xi[2])),
                      operator()(V3i(xi[0], 2*rad-1-xi[1], xi[2])),
                      operator()(V3i(xi[0], xi[1], 2*rad-1-xi[2])));

    }

    v3.resize(xis.size());
    for (size_t i = 0; i < xis.size(); ++i) {
        v3[i] = xis[i];
    }
}

template <>
inline XiMesh<Cylindrical>::XiMesh(int rad, double cut, const Vm v_ ) : 
        rad(rad), cut_(cut), a_(cut/rad), vol_(2*M_PI*sqr(a_)), shift(vm2vd(v_))
{
    xyzmap.resize(2*sqr(rad));
    int i = 0;
    for (int i1 = 0; i1 < 2*rad; ++i1)
        for (int i2 = 0; i2 < rad; ++i2) {
            V2i xi(i1, i2);
            V2d v = i2xi(xi) - shift;
            int j = flatten(xi);
            if (sqr(v) < sqr(cut_)) {
                xis.push_back(v + shift);
                vis.push_back(xi);
                xyzmap[j] = i;
                ++i;
            }
            else {
                xyzmap[j] = -1;
            }
        }

    mirr.resize(vis.size());
    for (size_t i = 0; i < vis.size(); ++i) {
        V2i xi(vis[i]);
        mirr[i] = operator()(V2i(2*rad-1-xi[0], xi[1]));
    }

    v3.resize(xis.size());
    for (size_t i = 0; i < xis.size(); ++i) {
        v3[i] = xis[i];
    }
}

#endif
