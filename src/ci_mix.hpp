#ifndef _CI_MIX_HPP_
#define _CI_MIX_HPP_

#include <vector>

#include "ci.hpp"
#include "ci_mixture.hpp"
#include "ci_multi.hpp"

#include "ximesh_mix.hpp"

#include "stencil.hpp"

enum TimeScheme {
    Continues,
    Euler
};
template <TimeScheme> struct trait_timescheme {};

inline std::ostream& operator<<(std::ostream& to, TimeScheme timescheme) {
    if       (timescheme == Continues) return to << "Continues";
    else if  (timescheme == Euler)     return to << "Euler"    ;
    else                               return to << timescheme ;
}

template<Symmetry symmetry>
std::tuple<V3d, double, V3d, V3d, V3d, V3d>
calcNodeCollide(const int i1, const int i2,
                const V3d v, const V3d w,
                const V3d n,
                const XiMeshMix<symmetry>& ximesh)
{
    return calcNodeCollideMixture(i1, i2, v, w, n, ximesh);
}

inline std::tuple<XiMeshMix<Cartesian>::Vi,
                  XiMeshMix<Cartesian>::Vi,
                  XiMeshMix<Cartesian>::Vd,
                  XiMeshMix<Cartesian>::Vd>
        calcNodeAfter(const int i1, const int i2,
                      const V3d v2, const V3d w2, 
                      const XiMeshMix<Cartesian>& ximesh)
{
    const XiMeshMix<Cartesian>::Vi vj = 
            ximesh.xi2i(XiMeshMix<Cartesian>::Vx(v2, ximesh.i2ci(i1)));
    const XiMeshMix<Cartesian>::Vi wj =
            ximesh.xi2i(XiMeshMix<Cartesian>::Vx(w2, ximesh.i2ci(i2)));
    return std::make_tuple(vj, wj, v2, w2);
}

inline std::tuple<XiMeshMix<Cylindrical>::Vi,
                  XiMeshMix<Cylindrical>::Vi,
                  XiMeshMix<Cylindrical>::Vd,
                  XiMeshMix<Cylindrical>::Vd>
        calcNodeAfter(const int i1, const int i2,
                      const V3d v2, const V3d w2, 
                      const XiMeshMix<Cylindrical>& ximesh)
{
    const V2d v1(v2[0], norm(V2d(v2[1], v2[2])));
    const V2d w1(w2[0], norm(V2d(w2[1], w2[2])));
    const XiMeshMix<Cylindrical>::Vi vj = 
            ximesh.xi2i(XiMeshMix<Cylindrical>::Vx(v1, ximesh.i2ci(i1)) );
    const XiMeshMix<Cylindrical>::Vi wj = 
            ximesh.xi2i(XiMeshMix<Cylindrical>::Vx(w1, ximesh.i2ci(i2)) );
    return std::make_tuple(vj, wj, v1, w1);
}

template <Symmetry symmetry>
double normRSwarm(const int i, const XiMeshMix<symmetry>& ximesh)
{
    return ximesh.a(i);
}

template <PowMethod powmethod, typename F, typename LogF>
void findLogf(const F& f, LogF& logf) {
    for (size_t i = 0; i < f.size(); i += 2) 
        log_<powmethod>(&f[i], &logf[i]);
    for (size_t i = 0; i < f.size(); ++i)
        if (f[i] < 0)
            std::cout << "f[" << i << "] = " << f[i] << std::endl;
        else if (f[i] == 0)
            logf[i] = std::numeric_limits<double>::min_exponent;
}

template <PowMethod powmethod, typename DF, class Nodes>
void ciIterMixEuler(DF& f, const Nodes& nodes)
{
    static std::vector<double> logf(f.size());
    findLogf<powmethod>(f, logf);
    static std::vector<double> fold(f.size());
    copy(f, fold);

    for (typename Nodes::const_iterator p = nodes.begin(); p != nodes.end(); ++p) {
        typename Nodes::const_reference n = *p;
        typedef typename Nodes::value_type Node;
        typedef typename Node::Sd Sd;

        const double fi1 = fold[n.i1];
        const double fi2 = fold[n.i2];

        const Sd logf1 = fAt(logf, n.j1);
        const Sd logf2 = fAt(logf, n.j2);

        const Sd r1 = n.x1; 
        const Sd r2 = n.x2; 

        const double g1 = interpF(logf1, r1);
        const double g2 = interpF(logf2, r2);

        const double ff1 = aff<Node::symmetry>(exp<powmethod>(g1, g2), n);

        const double ff = fi1 * fi2;

        const double d = ( - ff1 + ff ) * n.c; 	

        const Sd d1 = r1 * d;
        const Sd d2 = r2 * d;

        const Sd fbackup1 = fAt(f, n.j1);
        const Sd fbackup2 = fAt(f, n.j2);
        const double fibackup1 = f[n.i1];
        const double fibackup2 = f[n.i2];

        addF(f, n.j1, d1);
        addF(f, n.j2, d2);
        f[n.i1] -= d;
        f[n.i2] -= d;
        
        if ( lessZero(fAt(f, n.j1)) ||
             lessZero(fAt(f, n.j2)) ||
            (f[n.i1] < 0)           ||
            (f[n.i2] < 0)           )
        {
            fIs(f, n.j1, fbackup1);
            fIs(f, n.j2, fbackup2);
            f[n.i1] = fibackup1;
            f[n.i2] = fibackup2;
        }
    }
}

template <Symmetry symmetry, Volume volume = Symmetric, TimeScheme timescheme = Continues>
class ColliderMix {
    public:
        typedef XiMeshMix<symmetry> XiMeshMixType;
        typedef XiMesh<symmetry> XiMeshType;

        ColliderMix(const XiMeshMixType& ximesh) : ximesh_(ximesh) {
            std::cout << "Mix: symmetry = "       << symmetry     << 
                              ", volume = "       << volume       <<
                              ",\ntimescheme = "  << timescheme   << 
                              ", sizeof(Node) = " << sizeof(Node) << std::endl;
        }

        void gen(const double time_step, const int p,
                 const XiMeshMixType& ximesh,
                 const SimpleSection* section)
        {
            ciGen(time_step, p, ximesh, nc, section);
        }

        template <typename DF>
        double iter(DF& f) const {
            const PowMethod powmethod = FastSSE;
            if (timescheme == Continues)
                return ciIterMultiCont<powmethod>(f, nc, ximesh_);
            if (timescheme == Euler)
                ciIterMixEuler<powmethod>(f, nc);
            return 0;
        }

    private:
        typedef CollisionNodeMulti<symmetry, volume> Node;

        typedef typename XiMeshMixType::Vi Vi;
        typedef typename XiMeshMixType::Vd Vd;
        typedef Vd Vx;

        std::vector<Node> nc;
        const XiMeshMixType& ximesh_;
        
};

#endif
