#ifndef _CI_MIXTURE_HPP_
#define _CI_MIXTURE_HPP_

#include <vector>

#include "ci.hpp"
#include "ximesh_mixture.hpp"

template <Symmetry symmetry>
class ColliderMixture {
    public:
        typedef XiMeshMixture<symmetry> XiMeshMixtureType;
        typedef XiMesh<symmetry> XiMeshType;

        ColliderMixture(const XiMeshMixtureType& ximesh) : ximesh_(ximesh) {
            std::cout << "Mixture: symmetry = " << symmetry << std::endl;
        }

        void gen(const double time_step, const int p,
                 const XiMeshMixtureType& ximesh,
                 const SimpleSection* section) {
            ciGen(time_step, p, ximesh, nc, section);
        }

        template <typename DF>
        double iter(DF& f) const {
            const Interpolation interp = PowerInterp;
            return CiIter<interp, DF, NodeContainer>()(f, nc);
        }

    private:
        typedef CollisionNode<symmetry> Node;

        typedef typename XiMeshMixtureType::Vi Vi;
        typedef typename XiMeshMixtureType::Vd Vd;
        typedef Vd Vx;
        typedef std::vector<Node> NodeContainer;

        NodeContainer nc;
        const XiMeshMixtureType& ximesh_;
};

inline std::tuple<bool, double, V3d, V3d, V3i>
        calcNodeAfter(const int i1, const int i2,
                      const V3d v2, const V3d w2, 
                      const V3d o,  const V3d nn,
                      const XiMeshMixture<Cartesian>& ximesh)
{
    const XiMeshMixture<Cartesian>::Vi vi = 
            ximesh.xi2i(XiMeshMixture<Cartesian>::Vx(v2, ximesh.i2ci(i1)));
    const XiMeshMixture<Cartesian>::Vi wi =
            ximesh.xi2i(XiMeshMixture<Cartesian>::Vx(w2, ximesh.i2ci(i2)));

    if ((ximesh(vi) < 0) || (ximesh(wi) < 0))
        return std::make_tuple(false, 0.0, V3d(0.0), V3d(0.0), V3i(0));

    const double E = sqr(nn);
    const V3i xi = wi;

    return std::make_tuple(true, E, o * ximesh.m(i2), nn, xi);
}

inline std::tuple<bool, double, V2d, V2d, V3i>
        calcNodeAfter(const int i1, const int i2,
                      const V3d v2, const V3d w2, 
                      const V3d o,  const V3d nn,
                      const XiMeshMixture<Cylindrical>& ximesh)
{
    const V2d v1(v2[0], norm(V2d(v2[1], v2[2])));
    const V2d w1(w2[0], norm(V2d(w2[1], w2[2])));

    const XiMeshMixture<Cylindrical>::Vi vi = 
            ximesh.xi2i(XiMeshMixture<Cylindrical>::Vx(v1, ximesh.i2ci(i1)) );
    const XiMeshMixture<Cylindrical>::Vi wi = 
            ximesh.xi2i(XiMeshMixture<Cylindrical>::Vx(w1, ximesh.i2ci(i2)) );

    if ((ximesh(vi) < 0) || (ximesh(wi) < 0))
        return std::make_tuple(false, 0.0, V2d(0.0), V2d(0.0), V3i(0));

    const V3i xi = V3d(wi.vi, vi.vi[1]);
    const double E = ximesh.e(i1) + ximesh.e(i2);

    return std::make_tuple(true, E, v1, w1, xi);
}


inline std::tuple<bool, 
                  XiMeshMixture<Cartesian>::Vi, 
                  XiMeshMixture<Cartesian>::Vi, 
                  double, double>
        fit(const V3i xi, 
            const XiMeshMixture<Cartesian>::Vi vi, 
            const XiMeshMixture<Cartesian>::Vi wi,
            const V3d x, const V3d y,
            const XiMeshMixture<Cartesian>& ximesh)
{
    XiMeshMixture<Cartesian>::Vi xi2(xi, wi.i);
    XiMeshMixture<Cartesian>::Vi xi1(vi.vi + wi.vi - xi, vi.i);

    if ((ximesh(xi1) < 0) || (ximesh(xi2) < 0))
        return std::make_tuple(false, xi1, xi2, 0.0, 0.0);

    V3d    n = ximesh.i2xi(xi2) - x;
    double e = sqr(n);
    double q = sqr(n - y);
    return std::make_tuple(true, xi1, xi2, e, q);
}

inline std::tuple<bool, 
                  XiMeshMixture<Cylindrical>::Vi, 
                  XiMeshMixture<Cylindrical>::Vi, 
                  double, double>
        fit(const V3i xi, 
            const XiMeshMixture<Cylindrical>::Vi vi, 
            const XiMeshMixture<Cylindrical>::Vi wi,
            const V2d x, const V2d y,
            const XiMeshMixture<Cylindrical>& ximesh)
{
    XiMeshMixture<Cylindrical>::Vi 
            xi2( V2i(xi[0], xi[1]), wi.i );
    XiMeshMixture<Cylindrical>::Vi 
            xi1( V2i(vi.vi[0] + wi.vi[0] - xi[0], xi[2]), vi.i );
    
    const int i1 = ximesh(xi1);
    const int i2 = ximesh(xi2);
    if ((i1 < 0) || (i2 < 0))
        return std::make_tuple(false, xi1, xi2, 0.0, 0.0);

    const double e = ximesh.e(i1) + ximesh.e(i2);

    V2d v(ximesh.i2xi(xi1));
    V2d w(ximesh.i2xi(xi2));

    const double q = sqr(v - x) + sqr(w - y);
    return std::make_tuple(true, xi1, xi2, e, q);
}

template<Symmetry symmetry, template <Symmetry> class XiMeshMixtureType>
std::tuple<V3d, double, V3d, V3d, V3d, V3d>
calcNodeCollideMixture(const int i1, const int i2,
                       const V3d v, const V3d w,
                       const V3d n,
                       const XiMeshMixtureType<symmetry>& ximesh)
{
    const double m1 = ximesh.m(i1);
    const double m2 = ximesh.m(i2);

    const V3d    u = w / m2 - v / m1;
    const double g = norm(u);
    const V3d   nn = m1 * m2 / (m1 + m2) * n * norm(u);
    const V3d    o = (v + w) / (m1 + m2);
    const V3d   v2 = m1 * o - nn;
    const V3d   w2 = m2 * o + nn;
    return std::make_tuple(u, g, o, nn, v2, w2);
}

template<Symmetry symmetry>
std::tuple<V3d, double, V3d, V3d, V3d, V3d>
calcNodeCollide(const int i1, const int i2,
                const V3d v, const V3d w,
                const V3d n,
                const XiMeshMixture<symmetry>& ximesh)
{
    return calcNodeCollideMixture(i1, i2, v, w, n, ximesh);
}


#endif
