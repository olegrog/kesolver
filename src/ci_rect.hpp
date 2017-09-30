#ifndef _Ci_RECT_HPP_
#define _Ci_RECT_HPP_

#include "ximesh_rect.hpp"

#include "ci.hpp"
#include "ci_multi.hpp"

template <Symmetry symmetry, Volume volume = Symmetric>
class ColliderRect {
    public:
        typedef XiMeshRect<symmetry> XiMeshType;

        ColliderRect(const XiMeshType& ximesh) : ximesh_(ximesh) {
            std::cout << "Rect: symmetry = "      << symmetry     << 
                             ", volume = "        << volume       <<
                             ", sizeof(Node) = "  << sizeof(Node) << std::endl;
        }

        void gen(const double time_step, const int p,
                 const XiMeshType& ximesh,
                 const SimpleSection* section) {
            ciGen(time_step, p, ximesh, nc, section);
        }

        template <typename DF>
        double iter(DF& f) const {
            const PowMethod powmethod = FastSSE;
            return ciIterMultiCont<powmethod>(f, nc, ximesh_);
        }

    private:
        typedef CollisionNodeMulti<symmetry, volume> Node;

        typedef typename XiMeshType::Vi Vi;
        typedef typename XiMeshType::Vd Vd;
        typedef Vd Vx;

        std::vector<Node> nc;
        const XiMeshType& ximesh_;
};

inline std::tuple<XiMeshRect<Cartesian>::Vi,
                    XiMeshRect<Cartesian>::Vi,
                    XiMeshRect<Cartesian>::Vd,
                    XiMeshRect<Cartesian>::Vd>
        calcNodeAfter(const int i1, const int i2,
                      const V3d v2, const V3d w2, 
                      const XiMeshRect<Cartesian>& ximesh)
{
    const XiMeshRect<Cartesian>::Vi vj = ximesh.xi2i(v2);
    const XiMeshRect<Cartesian>::Vi wj = ximesh.xi2i(w2);
    return std::make_tuple(vj, wj, v2, w2);
}

inline std::tuple<XiMeshRect<Cylindrical>::Vi,
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
    return std::make_tuple(vj, wj, v1, w1);
}

template<Symmetry symmetry>
std::tuple<V3d, double, V3d, V3d, V3d, V3d>
calcNodeCollide(const int i1, const int i2,
                   const V3d v, const V3d w,
                   const V3d n,
                   const XiMeshRect<symmetry>& ximesh)
{
    return calcNodeCollideSimple(i1, i2, v, w, n, ximesh);
}

template <Symmetry symmetry>
double normRSwarm(const int i, const XiMeshRect<symmetry>& ximesh)
{
    return 1.;
}

#endif
