#ifndef _Ci_RECT_HPP_
#define _Ci_RECT_HPP_

#include "ximesh_rect.hpp"

#include "ci.hpp"
#include "ci_multi.hpp"

template <Symmetry symmetry, Volume volume = Symmetric>
class ColliderRect {
    public:
        typedef XiMeshRect<symmetry> XiMeshType;

        ColliderRect() {
            std::cout << "Rect: symmetry = "      << symmetry     << 
                             ", volume = "        << volume       <<
                             ", sizeof(Node) = "  << sizeof(Node) << std::endl;
        }

        void gen(const double time_step, const int p,
                 const XiMeshType& ximesh) {
            ciGen(time_step, p, ximesh, nc);
        }

        template <typename DF>
        void iter(DF& f) const {
            const PowMethod powmethod = FastSSE;
            ciIterMultiCont<powmethod>(f, nc);
        }

    private:
        typedef CollisionNodeMulti<symmetry, volume> Node;

        typedef typename XiMeshType::Vi Vi;
        typedef typename XiMeshType::Vd Vd;
        typedef Vd Vx;

        std::vector<Node> nc;
        
};

inline boost::tuple<XiMeshRect<Cartesian>::Vi,
                    XiMeshRect<Cartesian>::Vi,
                    XiMeshRect<Cartesian>::Vd,
                    XiMeshRect<Cartesian>::Vd>
        calcNodeAfter(const int i1, const int i2,
                      const V3d v2, const V3d w2, 
                      const XiMeshRect<Cartesian>& ximesh)
{
    const XiMeshRect<Cartesian>::Vi vj = ximesh.xi2i(v2);
    const XiMeshRect<Cartesian>::Vi wj = ximesh.xi2i(w2);
    return boost::make_tuple(vj, wj, v2, w2);
}

inline boost::tuple<XiMeshRect<Cylindrical>::Vi,
                    XiMeshRect<Cylindrical>::Vi,
                    XiMeshRect<Cylindrical>::Vd,
                    XiMeshRect<Cylindrical>::Vd>
        calcNodeAfter(const int i1, const int i2,
                      const V3d v2, const V3d w2, 
                      const XiMeshRect<Cylindrical>& ximesh)
{
    const V2d v1(v2[0], norm(V2d(v2[1], v2[2])));
    const V2d w1(w2[0], norm(V2d(w2[1], w2[2])));
    const XiMeshRect<Cylindrical>::Vi vj = ximesh.xi2i(v1);
    const XiMeshRect<Cylindrical>::Vi wj = ximesh.xi2i(w1);
    return boost::make_tuple(vj, wj, v1, w1);
}

template<Symmetry symmetry>
boost::tuple<V3d, double, V3d, V3d, V3d, V3d>
calcNodeCollide(const int i1, const int i2,
                   const V3d v, const V3d w,
                   const V3d n,
                   const XiMeshRect<symmetry>& ximesh)
{
    return calcNodeCollideSimple(i1, i2, v, w, n, ximesh);
}

inline
std::vector< XiMeshRect<Cartesian>::Vi >
stencilMulti(XiMeshRect<Cartesian>::Vi vj, const XiMeshRect<Cartesian>& )
{
    typedef XiMeshRect<Cartesian>::Vi Vi;
    std::vector<Vi> stencil;
    for (int s1 = -1; s1 < 2; ++s1)
        for (int s2 = -1; s2 < 2; ++s2)
            for (int s3 = -1; s3 < 2; ++s3) {
                const int s = std::abs(s1) + std::abs(s2) + std::abs(s3);
                if ( (s > 0) && (s < 3) ) 
                    stencil.push_back(vj + V3i(s1, s2, s3));
            }
    return stencil;
}

inline
std::vector< XiMeshRect<Cylindrical>::Vi >
stencilMulti(XiMeshRect<Cylindrical>::Vi vj, const XiMeshRect<Cylindrical>& )
{
    typedef XiMeshRect<Cylindrical>::Vi Vi;
    std::vector<Vi> stencil;
    for (int s1 = -1; s1 < 2; ++s1)
        for (int s2 = -1; s2 < 2; ++s2)
            if (std::abs(s1) + std::abs(s2) > 0)
                stencil.push_back(vj + V2i(s1, s2));
    return stencil;
}

template <Volume volume, Symmetry symmetry> 
Stencil<symmetry, volume, double> makeRSwarm(const typename SymmetryTrait<symmetry>::Vd x,
                                             const Stencil<symmetry, volume, int> j,
                                             const XiMeshRect<symmetry> mesh)
{
    return doMakeRSwarm<volume>(x, j, mesh);
}

template <Symmetry symmetry>
double normRSwarm(const int i, const XiMeshRect<symmetry> ximesh)
{
    return 1.;
}

#endif