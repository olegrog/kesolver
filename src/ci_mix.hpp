#ifndef _CI_MIX_HPP_
#define _CI_MIX_HPP_

#include <vector>

#include "ci.hpp"
#include "ci_mixture.hpp"
#include "ximesh_mix.hpp"
#include "stencil.hpp"

template <Symmetry symmetry, Volume volume>
struct XorR {
    typedef typename SymmetryTrait<symmetry>::Vd Vd;
    typedef Stencil<symmetry, volume, double> Xd;
    typedef Stencil<symmetry, volume, double> Sd;

    const Xd getX(const Vd x) const {
        return makeRSwarm<volume>(x);
    }
    const Sd getR(const Xd x) const {
        return x;
    }
};

template <Symmetry symmetry>
struct XorR<symmetry, Wide> {
    typedef typename SymmetryTrait<symmetry>::Vd Vd;
    typedef typename SymmetryTrait<symmetry>::Vd Xd;
    typedef Stencil<symmetry, Wide, double> Sd;

    const Xd getX(const Vd x) const {
        return x;
    }
    const Sd getR(const Xd x) const {
        return makeRSwarm<Wide>(x);
    }
};

template <Volume volume, Symmetry symmetry>
const typename XorR<symmetry, volume>::Xd getX(const typename XorR<symmetry, volume>::Vd x) {
    return XorR<symmetry, volume>().getX(x);
}

template <Volume volume, Symmetry symmetry>
const typename XorR<symmetry, volume>::Sd getR(const typename XorR<symmetry, volume>::Xd x) {
    return XorR<symmetry, volume>().getR(x);
}

template <Symmetry symmetry, Volume volume>
struct CollisionNodeMix;

template <Volume volume_>
struct CollisionNodeMix<Cartesian, volume_> {
    static const Volume   volume   = volume_;
    static const Symmetry symmetry = Cartesian;
    typedef Stencil<symmetry, volume_, int>    Si;
    typedef Stencil<symmetry, volume_, double> Sd;
    typedef typename XorR<symmetry, volume_>::Xd Xd;
    int i1, i2;
    Si  j1, j2;
    Xd  x1, x2;
    double c;
};

template <Volume volume_>
struct CollisionNodeMix<Cylindrical, volume_> {
    static const Volume   volume   = volume_;
    static const Symmetry symmetry = Cylindrical;
    typedef Stencil<symmetry, volume_, int>    Si;
    typedef Stencil<symmetry, volume_, double> Sd;
    typedef typename XorR<symmetry, volume_>::Xd Xd;
    int i1, i2;
    Si  j1, j2;
    Xd  x1, x2;
    double c, a;
};

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

inline boost::tuple<XiMeshMix<Cartesian>::Vi,
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
    return boost::make_tuple(vj, wj, v2, w2);
}

inline boost::tuple<XiMeshMix<Cylindrical>::Vi,
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
    return boost::make_tuple(vj, wj, v1, w1);
}

template <Volume volume>
const CollisionNodeMix<Cartesian, volume>
makeNode(const int i1 , const int i2,
         const Stencil<Cartesian, volume, int> j1,
         const Stencil<Cartesian, volume, int> j2,
         const V3d x1, const V3d x2, 
         const XiMeshMix<Cartesian>& ximesh)
{
    CollisionNodeMix<Cartesian, volume> node;
    node.x1 = getX<volume, Cartesian>(x1);
    node.x2 = getX<volume, Cartesian>(x2);
    return node;
}

template <Volume volume>
const CollisionNodeMix<Cylindrical, volume>
makeNode(const int i1 , const int i2,
         const Stencil<Cylindrical, volume, int> j1,
         const Stencil<Cylindrical, volume, int> j2,
         const V2d x1, const V2d x2, 
         const XiMeshMix<Cylindrical>& ximesh)
{
    CollisionNodeMix<Cylindrical, volume> node;
    Stencil<Cylindrical, volume, double> r1 = makeRSwarm<volume>(x1);
    Stencil<Cylindrical, volume, double> r2 = makeRSwarm<volume>(x2);
    node.a  = ximesh.vol(i1) * ximesh.vol(i2) / 
              (interpVol(j1, r1, ximesh) * interpVol(j2, r2, ximesh));
    node.x1 = getX<volume, Cylindrical>(x1);
    node.x2 = getX<volume, Cylindrical>(x2);
    return node;
}

inline
std::vector< XiMeshMix<Cartesian>::Vi >
stencilMix(XiMeshMix<Cartesian>::Vi vj)
{
    typedef XiMeshMix<Cartesian>::Vi Vi;
    std::vector<Vi> stencil;
    for (int s1 = -1; s1 < 2; ++s1)
        for (int s2 = -1; s2 < 2; ++s2)
            for (int s3 = -1; s3 < 2; ++s3) {
                const int s = std::abs(s1) + std::abs(s2) + std::abs(s3);
                if ( (s > 0) && (s < 3) ) 
                    stencil.push_back(Vi(vj.vi + V3i(s1, s2, s3), vj.i));
            }
    return stencil;
}

inline
std::vector< XiMeshMix<Cylindrical>::Vi >
stencilMix(XiMeshMix<Cylindrical>::Vi vj)
{
    typedef XiMeshMix<Cylindrical>::Vi Vi;
    std::vector<Vi> stencil;
    for (int s1 = -1; s1 < 2; ++s1)
        for (int s2 = -1; s2 < 2; ++s2)
            if (std::abs(s1) + std::abs(s2) > 0)
                stencil.push_back(Vi(vj.vi + V2i(s1, s2), vj.i));
    return stencil;
}


template <Volume volume, Symmetry symmetry>
boost::tuple< bool, typename CollisionNodeMix<symmetry, volume>::Si, 
                  typename XiMeshMix<symmetry>::Vd >
searchMix(typename XiMeshMix<symmetry>::Vd v1, typename XiMeshMix<symmetry>::Vi vj,
          const XiMeshMix<symmetry>& ximesh)
{
    typedef typename XiMeshMix<symmetry>::Vi Vi;
    typedef typename XiMeshMix<symmetry>::Vd Vd;
    typedef typename CollisionNodeMix<symmetry, volume>::Si Si;
    Vd x1 = (v1 - ximesh.i2xi(vj).v);
    Si j1 = makeISwarm<volume>(vj, x1, ximesh);
    if (lessZero(j1)) {
        bool b = false;
        double q = std::numeric_limits<double>::max();
        const std::vector<Vi> stencil = stencilMix(vj);
        for (typename std::vector<Vi>::const_iterator p = stencil.begin(); p != stencil.end(); ++p)
        {
            const Vi& vj1 = *p;
            const Vd x = (v1 - ximesh.i2xi(vj1).v);
            const double q1 = sqr(x);
            if (q1 < q) {
                Si j = makeISwarm<volume>(vj1, x, ximesh);
                if ( not lessZero(j) ) {
                    q = q1;
                    x1 = x;
                    j1 = j;
                    b = true;
                }
            }
        }
        return boost::make_tuple(b, j1, x1);
    }
    return boost::make_tuple(true, j1, x1);
}


template <Volume volume, Symmetry symmetry>
CollisionType calcNode(const korobov::Point& p, 
                       const XiMeshMix<symmetry>& ximesh,
                       CollisionNodeMix<symmetry, volume>& node,
                       const SimpleSection* section)
{
    typedef typename XiMeshMix<symmetry>::Vi Vi;
    typedef typename XiMeshMix<symmetry>::Vd Vd;
    typedef CollisionNodeMix<symmetry, volume> Node;

    int i1, i2;
    Vi vi, wi;
    V3d v, w, n;
    boost::tie(i1, i2, vi, wi, v, w, n) = calcNodeBefore(p, ximesh);

    double g;
    V3d    u, o, nn, v2, w2;
    boost::tie(u, g, o, nn, v2, w2) = calcNodeCollide(i1, i2, v, w, n, ximesh);

    Vi vj, wj;
    Vd v1, w1;
    boost::tie(vj, wj, v1, w1) = calcNodeAfter(i1, i2, v2, w2, ximesh);

    int i3 = ximesh(vj);
    int i4 = ximesh(wj);
    if ((i3 < 0) || (i4 < 0))
        return BadCollision;

    typename Node::Si j1, j2;
    Vd x1, x2;
    bool b;
    boost::tie(b, j1, x1) = searchMix<volume>(v1, vj, ximesh);
    if (not b) return BadSearch;
    boost::tie(b, j2, x2) = searchMix<volume>(w1, wj, ximesh);
    if (not b) return BadSearch;

    x1 /= ximesh.a(i3);
    x2 /= ximesh.a(i4);

    node = makeNode(i1, i2, j1, j2, x1, x2, ximesh); 

    node.i1  = i1;
    node.i2  = i2;

    node.j1  = j1;
    node.j2  = j2;

    double m1 = ximesh.m(i1);
    double m2 = ximesh.m(i2);
    node.c = g * section->section(g * std::sqrt(m1 * m2 / (m1 + m2)),
                                  dot(u, n) / g,
                                  ximesh.i2ci(i1), ximesh.i2ci(i2));

    return GoodOne;
}

template <typename CollisionNodeType>
inline double aff(const double ff, const CollisionNodeType n, SymmetryTrait<Cartesian>) {
    return ff;
}

template <typename CollisionNodeType>
inline double aff(const double ff, const CollisionNodeType n, SymmetryTrait<Cylindrical>) {
    return n.a * ff;
}

template <Symmetry symmetry, typename CollisionNodeType>
inline double aff(const double ff, const CollisionNodeType n)
{
    return aff(ff, n, SymmetryTrait<symmetry>());
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

template <PowMethod powmethod> inline
double exp(const double g1, const double g2);

template <> inline
double exp<Std>(const double g1, const double g2) {
    return std::exp(g1 + g2);
}

template <> inline
double exp<FastSSE>(const double g1, const double g2) {
    sse::d2_t x;
    x.d[0] = g1;
    x.d[1] = g2;
    x = sse::exp(x);
    return x.d[0] * x.d[1];
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

        const Sd r1 = getR<Node::volume, Node::symmetry>(n.x1); 
        const Sd r2 = getR<Node::volume, Node::symmetry>(n.x2); 

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

template <PowMethod powmethod, typename DF, class Nodes>
void ciIterMixCont(DF& f, const Nodes& nodes)
{
    int i1 = 0, i2 = 0;
    for (typename Nodes::const_iterator p = nodes.begin(); p != nodes.end(); ++p) {
        typename Nodes::const_reference n = *p;
        typedef typename Nodes::value_type Node;
        typedef typename Node::Sd Sd;

        const Sd f1 = fAt(f, n.j1);
        const Sd f2 = fAt(f, n.j2);
        const double fi1 = f[n.i1];
        const double fi2 = f[n.i2];

        const Sd logf1 = mylog<powmethod>(f1);
        const Sd logf2 = mylog<powmethod>(f2);

        const Sd r1 = getR<Node::volume, Node::symmetry>(n.x1); 
        const Sd r2 = getR<Node::volume, Node::symmetry>(n.x2); 

        const double g1 = interpF(logf1, r1);
        const double g2 = interpF(logf2, r2);

        const double ff1 = aff<Node::symmetry>(exp<powmethod>(g1, g2), n);
        const double ff = fi1 * fi2;

        const double d = ( - ff1 + ff ) * n.c; 	

        const Sd d1 = r1 * d;
        const Sd d2 = r2 * d;

        addF(f, n.j1, d1);
        addF(f, n.j2, d2);
        f[n.i1] -= d;
        f[n.i2] -= d;
        
        if ( lessZero(fAt(f, n.j1)) ||
             lessZero(fAt(f, n.j2)) ||
            (f[n.i1] < 0)           ||
            (f[n.i2] < 0)           )
        {
            fIs(f, n.j1, f1);
            fIs(f, n.j2, f2);
            f[n.i1] = fi1;
            f[n.i2] = fi2;
            ++i2;
        }
        ++i1;
    }
//    std::cout << "i1, i2 = " << i1 << ' ' << i2 << ' ' << (i2 + 0.0) / i1 << std::endl;
}

template <Symmetry symmetry, Volume volume = Symmetric, TimeScheme timescheme = Continues>
class ColliderMix {
    public:
        typedef XiMeshMix<symmetry> XiMeshMixType;
        typedef XiMesh<symmetry> XiMeshType;

        ColliderMix() {
            std::cout << "Mix: symmetry = "       << symmetry     << 
                              ", volume = "       << volume       <<
                              ",\ntimescheme = "  << timescheme   << 
                              ", sizeof(Node) = " << sizeof(Node) << std::endl;
        }

        void gen(const double time_step, const int p,
                 const XiMeshMixType& ximesh,
                 const SimpleSection* section) {
            ciGen(time_step, p, ximesh, nc, section);
        }

        template <typename DF>
        void iter(DF& f) const {
            const PowMethod powmethod = FastSSE;
            if (timescheme == Continues)
                ciIterMixCont<powmethod>(f, nc);
            else if (timescheme == Euler)
                ciIterMixEuler<powmethod>(f, nc);
        }

    private:
        typedef CollisionNodeMix<symmetry, volume> Node;

        typedef typename XiMeshMixType::Vi Vi;
        typedef typename XiMeshMixType::Vd Vd;
        typedef Vd Vx;

        std::vector<Node> nc;
        
};

#endif
