#ifndef _XIMESH_RECT_HPP_
#define _XIMESH_RECT_HPP_

#include <vector>
#include <iostream>
#include <cmath>

#include "ximesh.hpp"

template <Symmetry symmetry>
class XiMeshRect {
    public:
        typedef typename SymmetryTrait<symmetry>::Vm Vm;
        typedef typename SymmetryTrait<symmetry>::Vd Vd;
        typedef typename SymmetryTrait<symmetry>::Vi Vi;
        typedef typename SymmetryTrait<symmetry>::Vj Vj;
        typedef Vd Vx;
        typedef typename SymmetryTrait<symmetry>::VVd VVd;

        XiMeshRect(const VVd& vv, const VVd& vvol);

        size_t size() const { return xis.size(); }

        const Vi radius() const { return rad; }
        double cut() const { return max(cut_); }
        double vol() const { return vol_; }
        double vol(const int i) const 
                { return vols[i]; }
        double volLog(const int i) const
                { return vol_logs[i]; }

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

        const Vx i2xi(const Vi i) const;
        const Vi xi2i(const Vx xi) const;

        int mirror(int i, const int a) const;
        const Vj mirror(int i) const
                { return mirr[i]; }

        double vvi(int i, int a) const 
                { return vv[a][i]; }

        double vvoli(int i, int a) const 
                { return vvol[a][i2vi(i)[a]]; }

    private:
        Vi rad;
        Vd cut_;
        double vol_;
        VVd vv, vvol;
        std::vector<Vd>  xis;
        std::vector<Vi>  vis;
        std::vector<double> vols, vol_logs;
        std::vector<int> xyzmap;
        std::vector<Vj>  mirr;
        std::vector<V3d> v3;

        int flatten(const Vi i) const;
};

template <>
inline const XiMeshRect<Cartesian>::Vm XiMeshRect<Cartesian>::p(const int i) const {
    return xis[i];
}

template <>
inline const XiMeshRect<Cylindrical>::Vm XiMeshRect<Cylindrical>::p(const int i) const {
    return xis[i][0];
}

template <>
inline int XiMeshRect<Cartesian>::flatten(const Vi i) const {
    return dot(i, Vi(1, 2*rad[0], 2*rad[1]*2*rad[0]));
};

template <>
inline int XiMeshRect<Cylindrical>::flatten(const Vi i) const {
    return dot(i, Vi(1, 2*rad[0]));
};

template <>
inline int XiMeshRect<Cartesian>::operator()(const Vi i) const {
    int res = ((i >= 0) && (i < 2*rad)) ? xyzmap[flatten(i)] : -1;
    return res;
}

template <>
inline int XiMeshRect<Cylindrical>::operator()(Vi i) const {
    if (i[1] < 0) i[1] = - i[1] - 1;
//    std::cout << "i, rad = " << i << ' ' << rad << std::endl;
    return ((i >= 0) && (i < Vi(2*rad[0], rad[1]))) ? xyzmap[flatten(i)] : -1;
}

template <>
inline const XiMeshRect<Cartesian>::Vx XiMeshRect<Cartesian>::i2xi(const Vi i) const {
    return Vx(vv[0][i[0]], vv[1][i[1]], vv[2][i[2]]);
}

template <>
inline const XiMeshRect<Cylindrical>::Vx XiMeshRect<Cylindrical>::i2xi(const Vi i) const {
    return Vx(vv[0][i[0]], vv[1][i[1]]);
}

template <Symmetry symmetry>
inline const typename XiMeshRect<symmetry>::Vi XiMeshRect<symmetry>::xi2i(const Vx xi) const {
    Vi i;
    for (int j = 0; j < SymmetryTrait<symmetry>::ximesh_dim; ++j) {
        double x = xi[j];
        if (x < vv[j][0] - 0.5 * vvol[j][0]) {
            i[j] = -1;
        }
        else if (x > vv[j][vv[j].size()-1] + 0.5 * vvol[j][vv[j].size()-1]) {
            i[j] = vv[j].size();
        }
        else {
            int l = 0, r = vv[j].size(), c = 0;
            while (l < r) {
                c = (l + r) / 2;
                double   v  = vv[j][c],       vol = vvol[j][c];
                double   v1 = v - 0.5 * vol,  v2  = v + 0.5 * vol;
                if ((x >= v1) && (x < v2))
                    break;
                else if (x < v1)
                    r = c;
                else if (x >= v2)
                    l = c + 1;
            }
            i[j] = c;
        }
    }
    return i;
}

template <>
inline int XiMeshRect<Cartesian>::mirror(const int i, const int a) const {
    return mirr[i][a];
}

template <>
inline int  XiMeshRect<Cylindrical>::mirror(const int i, const int a) const {
    return mirr[i];
}


template <>
inline XiMeshRect<Cartesian>::XiMeshRect(const VVd& vv, const VVd& vvol) : 
        vol_(1.), vv(vv), vvol(vvol)
{
    rad  = Vi(vv[0].size() / 2, vv[1].size() / 2, vv[2].size() / 2);
    size_t size = 8 * rad[0] * rad[1] * rad[2];

    for (int i = 0; i < 3; ++i) {
        double x = 0.;
        for (size_t j = 0; j < vv[i].size(); ++j) {
            double y = std::fabs(vv[i][j]) + 0.5 * vvol[i][j];
            if (x < y) x = y;
        }
        cut_[i] = x;
    }
    
    std::cout << "vol_ = " << vol_ << std::endl;
    xyzmap.resize(size);
    int i = 0;
    for (int i1 = 0; i1 < 2 * rad[0]; ++i1)
        for (int i2 = 0; i2 < 2 * rad[1]; ++i2)
            for (int i3 = 0; i3 < 2 * rad[2]; ++i3) {
                V3i xi(i1, i2, i3);
                V3d v = i2xi(xi);
                int j = flatten(xi);
                if (sqr(v / cut_) < 1) {
                    xis.push_back(v);
                    vis.push_back(xi);
                    double vol = vvol[0][i1] * vvol[1][i2] * vvol[2][i3];
                    vols.push_back(vol);
                    vol_logs.push_back(std::log(vol));
                    xyzmap[j] = i;
                    ++i;
                }
                else 
                    xyzmap[j] = -1;
            }
    std::cout << "ximesh.size() = " << xis.size() << std::endl;
    mirr.resize(xis.size());
    for (size_t i = 0; i < vis.size(); ++i) {
        V3i xi(vis[i]);
        mirr[i] = V3i(operator()(V3i(2*rad[0]-1-xi[0], xi[1], xi[2])),
                      operator()(V3i(xi[0], 2*rad[1]-1-xi[1], xi[2])),
                      operator()(V3i(xi[0], xi[1], 2*rad[2]-1-xi[2])));
    }

    v3.resize(xis.size());
    for (size_t i = 0; i < xis.size(); ++i) {
        v3[i] = xis[i];
    }
}

template <>
inline XiMeshRect<Cylindrical>::XiMeshRect(const VVd& vv, const VVd& vvol) : 
        vol_(2*M_PI), vv(vv), vvol(vvol)
{
    rad  = Vi(vv[0].size() / 2, vv[1].size());
    size_t size = 2 * rad[0] * rad[1];

    for (int i = 0; i < 2; ++i) {
        double x = 0.;
        for (size_t j = 0; j < vv[i].size(); ++j) {
            double y = std::fabs(vv[i][j]) + 0.5 * vvol[i][j];
            if (x < y) x = y;
        }
        cut_[i] = x;
    }
    
    std::cout << "vol_ = " << vol_ << std::endl;
    xyzmap.resize(size);
    int i = 0;
    for (int i1 = 0; i1 < 2*rad[0]; ++i1)
        for (int i2 = 0; i2 < rad[1]; ++i2) {
            V2i xi(i1, i2);
            V2d v = i2xi(xi);
            int j = flatten(xi);
            if (sqr(v / cut_) < 1) {
                xis.push_back(v);
                vis.push_back(xi);
                double vol = v[1] * vvol[0][i1] * vvol[1][i2];
                vols.push_back(vol);
                vol_logs.push_back(std::log(vol));
                xyzmap[j] = i;
                ++i;
            }
            else {
                xyzmap[j] = -1;
            }
        }
    std::cout << "ximesh.size() = " << xis.size() << std::endl;
    mirr.resize(xis.size());
    for (size_t i = 0; i < vis.size(); ++i) {
        V2i xi(vis[i]);
        mirr[i] = operator()(V2i(2*rad[0]-1-xi[0], xi[1]));
    }

    v3.resize(xis.size());
    for (size_t i = 0; i < xis.size(); ++i) {
        v3[i] = xis[i];
    }
}

#endif
