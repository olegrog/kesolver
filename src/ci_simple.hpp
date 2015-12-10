#ifndef _CI_SIMPLE_HPP_
#define _CI_SIMPLE_HPP_

#include "ci.hpp"


template <Interpolation interp, Symmetry symmetry>
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
            CiIter<interp, DF, NodeContainer>()(f, nc);
        }

    private:
        typedef CollisionNode<symmetry> Node;

        typedef typename SymmetryTrait<symmetry>::Vi Vi;
        typedef typename SymmetryTrait<symmetry>::Vd Vd;
        typedef Vd Vx;
        typedef std::vector<Node> NodeContainer;

        NodeContainer nc;

};

inline std::tr1::tuple<bool, double, V3d, V3d, V3i>
        calcNodeAfter(const int i1, const int i2,
                      const V3d v2, const V3d w2, 
                      const V3d o,  const V3d nn,
                      const XiMesh<Cartesian>& ximesh)
{
    const V3i vi = ximesh.xi2i(v2);
    const V3i wi = ximesh.xi2i(w2);
    if ((ximesh(vi) < 0) || (ximesh(wi) < 0))
        return std::tr1::make_tuple(false, 0.0, 0.0, 0.0, 0);

    const double E = sqr(nn);
    const V3i xi = wi;

    return std::tr1::make_tuple(true, E, o, nn, xi);
}

inline std::tr1::tuple<bool, double, V2d, V2d, V3i>
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
        return std::tr1::make_tuple(false, 0.0, 0.0, 0.0, 0);

    const V3i xi = V3d(wi, vi[1]);
    const double E = ximesh.e(i1) + ximesh.e(i2);

    return std::tr1::make_tuple(true, E, v1, w1, xi);
}

inline std::tr1::tuple<bool, V3i, V3i, double, double>
        fit(const V3i xi, const V3i vi, const V3i wi, const V3d x, const V3d y,
            const XiMesh<Cartesian>& ximesh)
{
    V3i xi2 = xi;
    V3i xi1 = vi + wi - xi;
    if ((ximesh(xi1) < 0) || (ximesh(xi2) < 0))
        return std::tr1::make_tuple(false, xi1, xi2, 0.0, 0.0);
    V3d n = ximesh.i2xi(xi2) - x;
    double e = sqr(n);
    double q = sqr(n - y);
    return std::tr1::make_tuple(true, xi1, xi2, e, q);
}

inline std::tr1::tuple<bool, V2i, V2i, double, double>
        fit(const V3i xi, const V2i vi, const V2i wi, const V2d x, const V2d y,
            const XiMesh<Cylindrical>& ximesh)
{
    V2i xi2(xi[0], xi[1]);
    V2i xi1(vi[0] + wi[0] - xi[0], xi[2]);
    
    const int i1 = ximesh(xi1);
    const int i2 = ximesh(xi2);
    if ((i1 < 0) || (i2 < 0))
        return std::tr1::make_tuple(false, xi1, xi2, 0.0, 0.0);

    const double e = ximesh.e(i1) + ximesh.e(i2);

    V2d v(ximesh.i2xi(xi1));
    V2d w(ximesh.i2xi(xi2));

    const double q = sqr(v - x) + sqr(w - y);
    return std::tr1::make_tuple(true, xi1, xi2, e, q);
}

template<Symmetry symmetry>
std::tr1::tuple<V3d, double, V3d, V3d, V3d, V3d>
calcNodeCollide(const int i1, const int i2,
                   const V3d v, const V3d w,
                   const V3d n,
                   const XiMesh<symmetry>& ximesh)
{
    return calcNodeCollideSimple(i1, i2, v, w, n, ximesh);
}

#endif
