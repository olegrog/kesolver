#ifndef _Ci_MULTI_HPP_
#define _Ci_MULTI_HPP_

#include "stencil.hpp"

template <Symmetry symmetry, Volume volume>
struct CollisionNodeMulti;

template <Volume volume_>
struct CollisionNodeMulti<Cartesian, volume_> {
    static const Volume   volume   = volume_;
    static const Symmetry symmetry = Cartesian;
    typedef Stencil<symmetry, volume_, int>    Si;
    typedef Stencil<symmetry, volume_, double> Sd;
    typedef Stencil<symmetry, volume_, double> Xd;
    int i1, i2;
    Si  j1, j2;
    Xd  x1, x2;
    double c;
};

template <Volume volume_>
struct CollisionNodeMulti<Cylindrical, volume_> {
    static const Volume   volume   = volume_;
    static const Symmetry symmetry = Cylindrical;
    typedef Stencil<symmetry, volume_, int>    Si;
    typedef Stencil<symmetry, volume_, double> Sd;
    typedef Stencil<symmetry, volume_, double> Xd;
    int i1, i2;
    Si  j1, j2;
    Xd  x1, x2;
    double c, a;
};

template <Volume volume, template <Symmetry symmetry> class XiMeshType>
const CollisionNodeMulti<Cartesian, volume>
makeNode(const int i1 , const int i2,
         const Stencil<Cartesian, volume, int> j1,
         const Stencil<Cartesian, volume, int> j2,
         const V3d x1, const V3d x2, 
         const XiMeshType<Cartesian>& ximesh)
{
    CollisionNodeMulti<Cartesian, volume> node;
    Stencil<Cartesian, volume, double> r1 = makeRSwarm<volume>(x1, j1, ximesh);
    Stencil<Cartesian, volume, double> r2 = makeRSwarm<volume>(x2, j2, ximesh);
    node.x1 = r1;
    node.x2 = r2;
    return node;
}

template <Volume volume, template <Symmetry symmetry> class XiMeshType>
const CollisionNodeMulti<Cylindrical, volume>
makeNode(const int i1 , const int i2,
         const Stencil<Cylindrical, volume, int> j1,
         const Stencil<Cylindrical, volume, int> j2,
         const V2d x1, const V2d x2, 
         const XiMeshType<Cylindrical>& ximesh)
{
    CollisionNodeMulti<Cylindrical, volume> node;
    Stencil<Cylindrical, volume, double> r1 = makeRSwarm<volume>(x1, j1, ximesh);
    Stencil<Cylindrical, volume, double> r2 = makeRSwarm<volume>(x2, j2, ximesh);
    node.a  = ximesh.vol(i1) * ximesh.vol(i2) / 
              (interpVol(j1, r1, ximesh) * interpVol(j2, r2, ximesh));
    node.x1 = r1;
    node.x2 = r2;
    return node;
}

template <Volume volume, Symmetry symmetry, template <Symmetry> class XiMeshType>
boost::tuple< bool, typename CollisionNodeMulti<symmetry, volume>::Si, 
                  typename XiMeshType<symmetry>::Vd >
searchMulti(typename XiMeshType<symmetry>::Vd v1, typename XiMeshType<symmetry>::Vi vj,
          const XiMeshType<symmetry>& ximesh)
{
    typedef typename XiMeshType<symmetry>::Vi Vi;
    typedef typename XiMeshType<symmetry>::Vd Vd;
    typedef typename CollisionNodeMulti<symmetry, volume>::Si Si;
    Vd x1 = (v1 - ximesh.i2xi(vj));
    Si j1 = makeISwarm<volume>(vj, x1, ximesh);
    if (lessZero(j1)) {
        bool b = false;
        double q = std::numeric_limits<double>::max();
        const std::vector<Vi> stencil = stencilMulti(vj, ximesh);
        for (typename std::vector<Vi>::const_iterator p = stencil.begin(); p != stencil.end(); ++p)
        {
            const Vi& vj1 = *p;
            const Vd x = (v1 - ximesh.i2xi(vj1));
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

template <typename Point, Volume volume, Symmetry symmetry, typename XiMeshType>
CollisionType calcNode(const Point& p, 
                       const XiMeshType& ximesh,
                       CollisionNodeMulti<symmetry, volume>& node,
                       const SimpleSection* section)
{
    typedef typename XiMeshType::Vi Vi;
    typedef typename XiMeshType::Vd Vd;
    typedef CollisionNodeMulti<symmetry, volume> Node;

    int i1, i2;
    Vi vi, wi;
    V3d v, w, n;
    boost::tie(i1, i2, vi, wi, v, w, n) = calcNodeBefore(p, ximesh);

//    std::cout << "v, w = " << v << ' ' << w << std::endl;

    double g;
    V3d    u, o, nn, v2, w2;
    boost::tie(u, g, o, nn, v2, w2) = calcNodeCollide(i1, i2, v, w, n, ximesh);

//    std::cout << "u, g = " << u << ' ' << g << std::endl;

    Vi vj, wj;
    Vd v1, w1;
    boost::tie(vj, wj, v1, w1) = calcNodeAfter(i1, i2, v2, w2, ximesh);

//    std::cout << "v1, w1 = " << v1 << ' ' << w1 << std::endl;

    int i3 = ximesh(vj);
    int i4 = ximesh(wj);

//    std::cout << "v2, w2 = " << v2 << ' ' << w2 << std::endl;
//    std::cout << "vj, wj = " << vj << ' ' << wj << std::endl;
//    std::cout << "i3, i4 = " << i3 << ' ' << i4 << std::endl;

    if ((i3 < 0) || (i4 < 0)) {
        return BadCollision;
    }

    typename Node::Si j1, j2;
    Vd x1, x2;
    bool b;
    boost::tie(b, j1, x1) = searchMulti<volume>(v1, vj, ximesh);
    if (not b) return BadSearch;
    boost::tie(b, j2, x2) = searchMulti<volume>(w1, wj, ximesh);
    if (not b) return BadSearch;

//    std::cout << "j1, j2 = " << j1 << ' ' << j2 << std::endl;
//    std::cout << "x1, x2 = " << x1 << ' ' << x2 << std::endl;

    x1 /= normRSwarm(i3, ximesh);
    x2 /= normRSwarm(i4, ximesh);

    node = makeNode(i1, i2, j1, j2, x1, x2, ximesh); 

    node.i1  = i1;
    node.i2  = i2;

    node.j1  = j1;
    node.j2  = j2;

    node.c = g * section->section(g, dot(u, n) / g);
    
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
void ciIterMultiCont(DF& f, const Nodes& nodes)
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

//        std::cout << "f1 = " << f1 << std::endl;
//        std::cout << "f2 = " << f2 << std::endl;

        const Sd r1 = n.x1; 
        const Sd r2 = n.x2; 

//        std::cout << "r1 = " << r1 << std::endl;
//        std::cout << "r2 = " << r2 << std::endl;

        const double g1 = interpF(logf1, r1);
        const double g2 = interpF(logf2, r2);

//        std::cout << "g1 = " << g1 << std::endl;
//        std::cout << "g2 = " << g2 << std::endl;
/*
        if ((g1 != g1) || (g2 != g2)) {
            LABEL

            std::cout << "f1 = " << f1 << std::endl;
            std::cout << "f2 = " << f2 << std::endl;

            std::cout << "r1 = " << r1 << std::endl;
            std::cout << "r2 = " << r2 << std::endl;

            exit(-1);
        }
*/
        const double ff1 = aff<Node::symmetry>(exp<powmethod>(g1, g2), n);
/*
        if (ff1 != ff1) {
            LABEL

            std::cout << "g1 = " << g1 << std::endl;
            std::cout << "g2 = " << g2 << std::endl;

            exit(-1);
        }
*/
        const double ff = fi1 * fi2;

        const double d = ( - ff1 + ff ) * n.c; 	
/*
        if (d != d) {
            LABEL

            std::cout << "ff1 = " << ff1 << std::endl;
            std::cout << "ff  = " << ff  << std::endl;

            exit(-1);
        }
*/
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

#endif
