#ifndef _CI_SIMPLE_HPP_
#define _CI_SIMPLE_HPP_

#include "ci.hpp"


template <Symmetry symmetry>
class ColliderSimple {
    public:
        typedef XiMesh<symmetry> XiMeshType;

        void gen(const double time_step, const int p,
                 const XiMeshType& ximesh,
                 const SimpleSection* section) {
            ciGen(time_step, p, ximesh, nc, section);
        }

        template <typename DF>
        void iter(DF& f) const {
            ciIter(f, nc);
        }

    private:
        typedef CollisionNode<symmetry> Node;

        typedef typename SymmetryTrait<symmetry>::Vi Vi;
        typedef typename SymmetryTrait<symmetry>::Vd Vd;
        typedef Vd Vx;

        std::vector<Node> nc;

};

inline boost::tuple<bool, double, V3d, V3d, V3i>
        calcNodeAfter(const int i1, const int i2,
                      const V3d v2, const V3d w2, 
                      const V3d o,  const V3d nn,
                      const XiMesh<Cartesian>& ximesh)
{
    const V3i vi = ximesh.xi2i(v2);
    const V3i wi = ximesh.xi2i(w2);
    if ((ximesh(vi) < 0) || (ximesh(wi) < 0))
        return boost::make_tuple(false, 0.0, 0.0, 0.0, 0);

    const double E = sqr(nn);
    const V3i xi = wi;

    return boost::make_tuple(true, E, o, nn, xi);
}

inline boost::tuple<bool, double, V2d, V2d, V3i>
        calcNodeAfter(const int i1, const int i2,
                      const V3d v2, const V3d w2, 
                      const V3d o,  const V3d nn,
                      const XiMesh<Cylindrical>& ximesh)
{
    const V2d v1(v2[0], norm(V2d(v2[1], v2[2])));
    const V2d w1(w2[0], norm(V2d(w2[1], w2[2])));

    const V2i vi = ximesh.xi2i(v1);
    const V2i wi = ximesh.xi2i(w1);
    if ((ximesh(vi) < 0) || (ximesh(wi) < 0))
        return boost::make_tuple(false, 0.0, 0.0, 0.0, 0);

    const V3i xi = V3d(wi, vi[1]);
    const double E = ximesh.e(i1) + ximesh.e(i2);

    return boost::make_tuple(true, E, v1, w1, xi);
}

inline boost::tuple<bool, V3i, V3i, double, double>
        fit(const V3i xi, const V3i vi, const V3i wi, const V3d x, const V3d y,
            const XiMesh<Cartesian>& ximesh)
{
    V3i xi2 = xi;
    V3i xi1 = vi + wi - xi;
    if ((ximesh(xi1) < 0) || (ximesh(xi2) < 0))
        return boost::make_tuple(false, xi1, xi2, 0.0, 0.0);
    V3d n = ximesh.i2xi(xi2) - x;
    double e = sqr(n);
    double q = sqr(n - y);
    return boost::make_tuple(true, xi1, xi2, e, q);
}

inline boost::tuple<bool, V2i, V2i, double, double>
        fit(const V3i xi, const V2i vi, const V2i wi, const V2d x, const V2d y,
            const XiMesh<Cylindrical>& ximesh)
{
    V2i xi2(xi[0], xi[1]);
    V2i xi1(vi[0] + wi[0] - xi[0], xi[2]);
    
    const int i1 = ximesh(xi1);
    const int i2 = ximesh(xi2);
    if ((i1 < 0) || (i2 < 0))
        return boost::make_tuple(false, xi1, xi2, 0.0, 0.0);

    const double e = ximesh.e(i1) + ximesh.e(i2);

    V2d v(ximesh.i2xi(xi1));
    V2d w(ximesh.i2xi(xi2));

    const double q = sqr(v - x) + sqr(w - y);
    return boost::make_tuple(true, xi1, xi2, e, q);
}

template<Symmetry symmetry>
boost::tuple<V3d, double, V3d, V3d, V3d, V3d>
calcNodeCollide(const int i1, const int i2,
                   const V3d v, const V3d w,
                   const V3d n,
                   const XiMesh<symmetry>& ximesh)
{
    const V3d    u = w - v;
    const double g = norm(u);
    const V3d    o = 0.5 * (v + w);
    const V3d   nn = 0.5 * n * g;
    const V3d   v2 = o - nn;
    const V3d   w2 = o + nn;
    return boost::make_tuple(u, g, o, nn, v2, w2);
}

#endif
