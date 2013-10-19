#ifndef _CI_HPP_
#define _CI_HPP_

#include <vector>
#include <algorithm>
#include <boost/tuple/tuple.hpp>
#include <cassert>

#include "symmetry.hpp"
#include "distribution_function.hpp"
#include "integr_grid.hpp"
#include "section.hpp"
#include "ximesh.hpp"
#include "random.hpp"
#include "sse.hpp"
#include "section.hpp"

template <Symmetry symmetry>
struct CollisionNode;

template <>
struct CollisionNode<Cartesian> {
    int i1, i2;
    int i1l, i2l;
    int i1m, i2m;
    int type;
    double c, r;
};

template <>
struct CollisionNode<Cylindrical> {
    int i1, i2;
    int i1l, i2l;
    int i1m, i2m;
    int type;
    double c, r, a;
};

enum CollisionType {
    BadCollision,
    BadSearch,
    SameXi,
    GoodOne,
    GoodTwo
};

template <Symmetry symmetry>
struct GridType;

template <>
struct GridType<Cartesian> {
    typedef korobov::Grid Grid;
};

template <>
struct GridType<Cylindrical> {
    typedef korobov::Grid Grid;
};

template <Symmetry symmetry, template <Symmetry symmetry> class XiMeshType, typename Nodes>
void ciGen(const double time_step, const int p,
           const XiMeshType<symmetry>& ximesh, Nodes& nc, 
           const SimpleSection* section)
{
    int ss[9];
    for (int i = 0; i < 9; ++i)
        ss[i] = 0;

    typedef typename GridType<symmetry>::Grid Grid;
    Grid grid(p);
    std::cout << grid.size() << std::endl;

    nc.clear();
    int n_nu = 0;
    for (typename Grid::iterator iter = grid.begin();
            iter != grid.end(); ++iter) {
        const typename Grid::point& point = *iter;

        CollisionType type;
        typename Nodes::value_type node;
        type = calcNode(point, ximesh, node, section);

        ss[type]++;
        if (type != BadSearch)
            ++n_nu;
        if ((type == GoodOne) || (type == GoodTwo))
            nc.push_back(node);
    }

    std::cout << "n_calc = " << nc.size() << " n_nu = " << n_nu << std::endl;
    for (int j = 0; j < 9; j++) 
        std::cout << ss[j] << ' ';
    std::cout << std::endl;

    const int nk = ximesh.size();
    const double B = (1/std::sqrt(2)/M_PI) * 4*M_PI * nk * nk * ximesh.vol() / n_nu / 4 * time_step;

    for (typename Nodes::iterator p = nc.begin(); p != nc.end(); ++p) {
        p->c *= B;
        double c = p->c;
        if (c != c) {
            LABEL
            std::cout << "c != c" << std::endl;
            exit(-1);
        }
    }

    std::random_shuffle(nc.begin(), nc.end());
}

template <typename Point, template <Symmetry symmetry> class XiMeshType>
inline boost::tuple<int, int, 
                  typename XiMeshType<Cartesian>::Vi,
                  typename XiMeshType<Cartesian>::Vi,
                  V3d, V3d,
                  V3d>
        calcNodeBefore(const Point& p,
                       const XiMeshType<Cartesian>& ximesh) 
{ 
    const int i1 = toInt(ximesh.size() * p[0]);
    const int i2 = toInt(ximesh.size() * p[1]);

    const typename XiMeshType<Cartesian>::Vi vi = ximesh.i2vi(i1);
    const typename XiMeshType<Cartesian>::Vi wi = ximesh.i2vi(i2);

    const V3d v = ximesh.i2xi(vi);
    const V3d w = ximesh.i2xi(wi);

    const V3d n = randomSphere(V2d(p[2], p[3]));

    return boost::make_tuple(i1, i2, vi, wi, v, w, n);
}

template <typename Point, template <Symmetry symmetry> class XiMeshType>
inline boost::tuple<int, int, 
                  typename XiMeshType<Cylindrical>::Vi,
                  typename XiMeshType<Cylindrical>::Vi,
                  V3d, V3d,
                  V3d>
        calcNodeBefore(const Point& p,
                       const XiMeshType<Cylindrical>& ximesh) 
{ 
    const int i1 = toInt(ximesh.size() * p[0]);
    const int i2 = toInt(ximesh.size() * p[1]);

    const typename XiMeshType<Cylindrical>::Vi vi = ximesh.i2vi(i1);
    const typename XiMeshType<Cylindrical>::Vi wi = ximesh.i2vi(i2);

    const V2d xi1 = ximesh.i2xi(vi);
    const V2d xi2 = ximesh.i2xi(wi);

    const double phi1 = 2 * M_PI * p[2];
    const double phi2 = 2 * M_PI * p[3];
    
    const V3d v  = rotate(V3d(xi1), V3d(1., 0., 0.), phi1);
    const V3d w  = rotate(V3d(xi2), V3d(1., 0., 0.), phi2);

    const V3d n = randomSphere(V2d(p[4], p[5]));

    return boost::make_tuple(i1, i2, vi, wi, v, w, n);
}

template<Symmetry symmetry, template <Symmetry symmetry> class XiMeshType>
boost::tuple<V3d, double, V3d, V3d, V3d, V3d>
calcNodeCollideSimple(const int i1, const int i2,
                      const V3d v, const V3d w,
                      const V3d n,
                      const XiMeshType<symmetry>& ximesh)
{
    const V3d    u = w - v;
    const double g = norm(u);
    const V3d    o = 0.5 * (v + w);
    const V3d   nn = 0.5 * n * g;
    const V3d   v2 = o - nn;
    const V3d   w2 = o + nn;
    return boost::make_tuple(u, g, o, nn, v2, w2);
}

template<Symmetry symmetry, template<Symmetry symmetry> class XiMeshType>
boost::tuple<int, double, 
           typename XiMeshType<symmetry>::Vi, typename XiMeshType<symmetry>::Vi,
           typename XiMeshType<symmetry>::Vi, typename XiMeshType<symmetry>::Vi
           >
        search(const std::vector<V3i>& stencil, const double E, 
               const typename XiMeshType<symmetry>::Vi vi, 
               const typename XiMeshType<symmetry>::Vi wi, 
               const typename XiMeshType<symmetry>::Vd x,
               const typename XiMeshType<symmetry>::Vd y,
               const XiMeshType<symmetry>& ximesh)
{
    typedef typename XiMeshType<symmetry>::Vi Vi;
    Vi xi2l, xi2m, xi1l, xi1m;
    double r = std::numeric_limits<double>::max();
    double q = std::numeric_limits<double>::max();
    for (std::vector<V3i>::const_iterator pl = stencil.begin(); pl != stencil.end(); ++pl) { 
        Vi xi2l_, xi1l_;
        double el, ql;
        bool bl;
        boost::tie(bl, xi1l_, xi2l_, el, ql) = fit(*pl, vi, wi, x, y, ximesh);
        if (bl) {
            if (el > E) {
                for (std::vector<V3i>::const_iterator pm = stencil.begin(); pm != stencil.end(); ++pm) {
                    Vi xi2m_, xi1m_;
                    double em, qm;
                    bool bm;
                    boost::tie(bm, xi1m_, xi2m_, em, qm) = fit(*pm, vi, wi, x, y, ximesh);
                    if (bm) {
                        if (em <= E+1e-12) {
                            double r_ = (E-el)/(em-el);
                            double q1 = (1-r_)*ql + r_*qm;
                            if (q1 < q) {
                                q = q1;
                                r = r_;
                                xi1l = xi1l_; xi1m = xi1m_;
                                xi2l = xi2l_; xi2m = xi2m_;
                                if (q < 1e-12)
                                    goto forend;
                            }
                        }
                    }
                }
            }
        }
    }
    forend:

    if (r == std::numeric_limits<double>::max()) 	
        return boost::make_tuple(0, 0.0, xi1l, xi2l, xi1m, xi2m);

    return boost::make_tuple(1, r, xi1l, xi2l, xi1m, xi2m);
}

template<Symmetry symmetry, template<Symmetry symmetry> class XiMeshType>
boost::tuple<int, double, 
           typename XiMeshType<symmetry>::Vi, typename XiMeshType<symmetry>::Vi,
           typename XiMeshType<symmetry>::Vi, typename XiMeshType<symmetry>::Vi
           >
        search2(const std::vector<V3i>& stencil, const double E, 
               const typename XiMeshType<symmetry>::Vi vi, 
               const typename XiMeshType<symmetry>::Vi wi, 
               const typename XiMeshType<symmetry>::Vd x,
               const typename XiMeshType<symmetry>::Vd y,
               const XiMeshType<symmetry>& ximesh)
{
    typedef typename XiMeshType<symmetry>::Vi Vi;
    Vi xi2l, xi2m, xi1l, xi1m;
    double r = std::numeric_limits<double>::max();
    double q = std::numeric_limits<double>::max();

    std::vector<V3i>::const_iterator pl = stencil.begin();

    Vi xi2l_, xi1l_;
    double el, ql;
    bool bl;
    boost::tie(bl, xi1l_, xi2l_, el, ql) = fit(*pl, vi, wi, x, y, ximesh);
    if (bl) {
        for (std::vector<V3i>::const_iterator pm = stencil.begin(); pm != stencil.end(); ++pm) {
            Vi xi2m_, xi1m_;
            double em, qm;
            bool bm;
            boost::tie(bm, xi1m_, xi2m_, em, qm) = fit(*pm, vi, wi, x, y, ximesh);
            if (bm) {
                double r_ = (E-el)/(em-el);
                if ( (r_ >= 0) && (r_ < 1) ) {
                    double q1 = (1-r_)*ql + r_*qm;
                    if (q1 < q) {
                        q = q1;
                        r = r_;
                        xi1l = xi1l_; xi1m = xi1m_;
                        xi2l = xi2l_; xi2m = xi2m_;
                        if (q < 1e-12)
                            break;
                    }
                }
            }
        }
    }

    if (r == std::numeric_limits<double>::max()) 	
        return boost::make_tuple(0, 0.0, xi1l, xi2l, xi1m, xi2m);

    return boost::make_tuple(1, r, xi1l, xi2l, xi1m, xi2m);
}

template<template<Symmetry symmetry> class XiMeshType>
inline const CollisionNode<Cartesian>
makeNode(const int iv , const int iw,
         const int ivl, const int iwl,
         const int ivm, const int iwm,
         const double r,
         const XiMeshType<Cartesian>& ximesh)
{
    CollisionNode<Cartesian> node;
    return node;
}

template<template<Symmetry symmetry> class XiMeshType>
inline const CollisionNode<Cylindrical>
makeNode(const int iv , const int iw,
         const int ivl, const int iwl,
         const int ivm, const int iwm,
         const double r,
         const XiMeshType<Cylindrical>& ximesh)
{
    CollisionNode<Cylindrical> node;
    const double x  = ximesh.vol(iv ) * ximesh.vol(iw );
    const double xl = ximesh.vol(ivl) * ximesh.vol(iwl);
    const double xm = ximesh.vol(ivm) * ximesh.vol(iwm);
    if (std::abs(1-r) > 1e-12)
        node.a = x / xl * pow(xl / xm, r);
    else
        node.a = x / xm;
    return node;
}

template <typename Point, Symmetry symmetry, template <Symmetry symmetry> class XiMeshType>
CollisionType calcNode(const Point& p, 
                       const XiMeshType<symmetry>& ximesh,
                       CollisionNode<symmetry>& node,
                       const SimpleSection* section)
{
    typedef typename XiMeshType<symmetry>::Vi Vi;
    typedef typename XiMeshType<symmetry>::Vd Vd;

    int i1, i2;
    Vi vi, wi;
    V3d v, w, n;
    boost::tie(i1, i2, vi, wi, v, w, n) = calcNodeBefore(p, ximesh);

    double g;
    V3d    u, o, nn, v2, w2;
    boost::tie(u, g, o, nn, v2, w2) = calcNodeCollide(i1, i2, v, w, n, ximesh);

    Vd x, y;
    double E;
    V3i xi;
    bool b;
    boost::tie(b, E, x, y, xi) = calcNodeAfter(i1, i2, v2, w2, o, nn, ximesh);

    if (!b)
        return BadCollision;

    std::vector<V3i> stencil;
    stencil.push_back(xi);
    for (int s1 = -1; s1 < 2; ++s1)
        for (int s2 = -1; s2 < 2; ++s2)
            for (int s3 = -1; s3 < 2; ++s3) {
                int s = abs(s1) + abs(s2) + abs(s3);
                if ( (s > 0) && (s <= 2) )
                    stencil.push_back(xi + V3i(s1, s2, s3));
            }

    double r;
    Vi xi1l, xi2l, xi1m, xi2m;
    int l;
    boost::tie(l, r, xi1l, xi2l, xi1m, xi2m) = 
            search2(stencil, E, vi, wi, x, y, ximesh);
    if (l == 0)
        return BadSearch;
 
    int i1l = ximesh( xi1l );
    int i1m = ximesh( xi1m );
    int i2l = ximesh( xi2l );
    int i2m = ximesh( xi2m );

    assert((i1l >= 0) && (i1m >= 0) && (i2l >= 0) && (i2m >= 0));
   
    if ( ( (i1 == i1l) && (i2 == i2l) ) || ( (i1 == i1m) && (i2 == i2m) ) )
        return SameXi;

    CollisionType type;

    node = makeNode(i1, i2, i1l, i2l, i1m, i2m, r, ximesh); 
    node.r = r;
    
    if (std::abs(1 - r) < 1e-12) {
        node.type = 1;
        type = GoodOne;
    }
    else if (std::abs(r) < 1e-12) {
        node.type = 1;
        i1m = i1l;
        i2m = i2l;
        type = GoodOne;
    }
    else {
        node.type = 2;
        type = GoodTwo;
    }

    node.i1  = i1;
    node.i1l = i1l;
    node.i1m = i1m;

    node.i2  = i2;
    node.i2l = i2l;
    node.i2m = i2m;

    node.c = g * section->section(g, dot(u, n) / g);

    return type;
}

inline double delta(const double fv, const double fw,
                    const double fv1, const double fw1, 
                    const CollisionNode<Cartesian> nc)
{
    return ( - fv1 * fw1 + fv * fw) * nc.c;
}

inline double delta(const double fv, const double fw,
                    const double fv1, const double fw1, 
                    const CollisionNode<Cylindrical> nc)
{
    return ( - nc.a * fv1 * fw1 + fv * fw) * nc.c;
}

template <typename DF, typename Nodes>
int ciIter(DF& f, const Nodes& nodes)
{
    int k = 0;
//    int l = 0;
    for (typename Nodes::const_iterator p = nodes.begin(); p != nodes.end(); ++p) {
        typename Nodes::const_reference n = *p;
        if (n.type == 2) {
            sse::d2_t x, y, z, w, v;

            x.d[0] = f[n.i1l];
            x.d[1] = f[n.i1m];

            z.d[0] = f[n.i2l];
            z.d[1] = f[n.i2m];

            w = sse::mul(x, z);
        
            y.d[0] = 1. - n.r;
            y.d[1] = n.r;
        
            v = sse::pow(w, y);

            double rr5 = f[n.i1];
            double rr6 = f[n.i2];
            double d = delta(rr5, rr6, v.d[0], v.d[1], *p);
            
            double dl = (1. - n.r) * d;
            double dm = n.r * d;

/*
            double d1 = delta(rr5, rr6, std::pow(f[n.i1l] * f[n.i2l], 1-n.r),
                                        std::pow(f[n.i1m] * f[n.i2m], n.r), *p);

            double h = d1 * ( std::log(f[n.i1])  + std::log(f[n.i2])   - 
                   (1-n.r) * (std::log(f[n.i1l]) + std::log(f[n.i2l])) -
                      n.r  * (std::log(f[n.i1m]) + std::log(f[n.i2m])) );

            if (h < 0) {
                std::cout << "h = " << h << std::endl;
            }
*/

            f[n.i1l] += dl;
            f[n.i2l] += dl;
            f[n.i1m] += dm;
            f[n.i2m] += dm;
            f[n.i1]  -= d;
            f[n.i2]  -= d;

            if ((f[n.i1l] < 0) ||
                (f[n.i1m] < 0) ||
                (f[n.i2l] < 0) ||
                (f[n.i2m] < 0) ||
                (f[n.i1]  < 0) || 
                (f[n.i2]  < 0))
            {
                f[n.i1l] = x.d[0];
                f[n.i1m] = x.d[1];
                f[n.i2l] = z.d[0];
                f[n.i2m] = z.d[1];
                f[n.i1 ] = rr5;
                f[n.i2 ] = rr6;
/*
                d = 2 * rr5 * rr6 * n.c;

                dl = (1. - n.r) * d;
                dm = n.r * d;

                f[n.i1l] += dl;
                f[n.i2l] += dl;
                f[n.i1m] += dm;
                f[n.i2m] += dm;
                f[n.i1]  -= d;
                f[n.i2]  -= d;

                if ((f[n.i1l] < 0) ||
                    (f[n.i1m] < 0) ||
                    (f[n.i2l] < 0) ||
                    (f[n.i2m] < 0) ||
                    (f[n.i1]  < 0) || 
                    (f[n.i2]  < 0))
                {
                    f[n.i1l] = x.d[0];
                    f[n.i1m] = x.d[1];
                    f[n.i2l] = z.d[0];
                    f[n.i2m] = z.d[1];
                    f[n.i1 ] = rr5;
                    f[n.i2 ] = rr6;

                    ++k;
                }
                else {
                    ++l;
                }
*/
            }
        }
        else {
            double g1 = f[n.i1];
            double g2 = f[n.i2];
            double g3 = f[n.i1m];
            double g4 = f[n.i2m];

            double d = delta(g1, g2, g3, g4, *p);

            f[n.i1]  -= d;
            f[n.i2]  -= d;
            f[n.i1m] += d;
            f[n.i2m] += d;

            if ((f[n.i1m] < 0) || 
                (f[n.i2m] < 0) || 
                (f[n.i1]  < 0) || 
                (f[n.i2]  < 0))
            {
                f[n.i1] = g1;
                f[n.i2] = g2;
                f[n.i1m] = g3;
                f[n.i2m] = g4;
                ++k;
            }
        }
    }
//    std::cout << "k / nodes.size() = " << static_cast<double>(k) / nodes.size() << std::endl;
//    std::cout << "l / nodes.size() = " << static_cast<double>(l) / nodes.size() << std::endl;
    return k;
}

#endif
