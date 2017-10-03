#ifndef _Ci_MULTI_HPP_
#define _Ci_MULTI_HPP_

#include "stencil.hpp"

template <Symmetry symmetry_, Volume volume_>
struct CollisionNodeMulti {
    static const Volume   volume   = volume_;
    static const Symmetry symmetry = symmetry_;
    typedef Stencil<symmetry, volume, int>    Si;
    typedef Stencil<symmetry, volume, double> Sd;
    typedef Stencil<symmetry, volume, double> Xd;
    int i1, i2;
    Si  j1, j2;
    Xd  x1, x2;     // projection weights
    Xd  w1, w2;     // interpolation weights
    double c, a;
};

template <Volume volume, Symmetry symmetry>
Stencil<symmetry, volume, double> createRSwarm(const typename SymmetryTrait<symmetry>::Vd x,
                                               const Stencil<symmetry, volume, int> j,
                                               const XiMeshMix<symmetry>& mesh)
{
    return makeRSwarm<volume>(x);
}

template <Volume volume, Symmetry symmetry>
Stencil<symmetry, volume, double> createWSwarm(const typename SymmetryTrait<symmetry>::Vd x,
                                               const Stencil<symmetry, volume, int> j,
                                               const XiMeshMix<symmetry>& mesh)
{
    return makeRSwarm<volume>(x);
}

template <Volume volume, Symmetry symmetry>
Stencil<symmetry, volume, double> createRSwarm(const typename SymmetryTrait<symmetry>::Vd x,
                                               const Stencil<symmetry, volume, int> j,
                                               const XiMeshRect<symmetry>& mesh)
{
    return doMakeRSwarm<volume>(x, j, mesh);
}

template <Volume volume, Symmetry symmetry, template <Symmetry> class XiMeshType>
struct MakeWSwarm;

template <Volume volume, template <Symmetry> class XiMeshType>
struct MakeWSwarm<volume, Cylindrical, XiMeshType> {
    static const Symmetry symmetry = Cylindrical;
    const Stencil<symmetry, volume, double> operator() (
        const typename SymmetryTrait<symmetry>::Vd x,
        const Stencil<symmetry, volume, int> j,
        const XiMeshType<symmetry>& ximesh) const
    {
        Stencil<symmetry, volume, double> s = {};
        double h1 = 0.5 * ( ximesh.vvoli(j.i, 0) + ximesh.vvoli(j.o, 0) );
        double h2 = 0.5 * ( ximesh.vvoli(j.j, 1) + ximesh.vvoli(j.o, 1) );
        V2d ax = abs(x);

        s.i = ax[0]/h1;
        s.j = ax[1]/h2;
        s.o = 1 - s.i - s.j;
        return s;
    }
};

template <Volume volume, template <Symmetry> class XiMeshType>
struct MakeWSwarm<volume, Cartesian, XiMeshType> {
    static const Symmetry symmetry = Cartesian;
    const Stencil<symmetry, volume, double> operator() (
        const typename SymmetryTrait<symmetry>::Vd x,
        const Stencil<symmetry, volume, int> j,
        const XiMeshType<symmetry>& ximesh) const
    {
        Stencil<symmetry, volume, double> s = {};
        double h1 = 0.5 * ( ximesh.vvoli(j.i, 0) + ximesh.vvoli(j.o, 0) );
        double h2 = 0.5 * ( ximesh.vvoli(j.j, 1) + ximesh.vvoli(j.o, 1) );
        double h3 = 0.5 * ( ximesh.vvoli(j.k, 2) + ximesh.vvoli(j.o, 2) );
        V3d ax = abs(x);

        s.i = ax[0]/h1;
        s.j = ax[1]/h2;
        s.k = ax[2]/h3;
        s.o = 1 - s.i - s.j - s.k;
        return s;
    }
};

template <Volume volume, Symmetry symmetry, template <Symmetry> class XiMeshType>
Stencil<symmetry, volume, double> createWSwarm(const typename SymmetryTrait<symmetry>::Vd x,
                                               const Stencil<symmetry, volume, int> j,
                                               const XiMeshType<symmetry>& ximesh)
{
    return MakeWSwarm<volume, symmetry, XiMeshType>()(x, j, ximesh);
}

template <Volume volume, template <Symmetry symmetry> class XiMeshType>
const CollisionNodeMulti<Cartesian, volume>
makeNode(const int i1 , const int i2,
         const Stencil<Cartesian, volume, int> j1,
         const Stencil<Cartesian, volume, int> j2,
         const V3d x1, const V3d x2, 
         const XiMeshType<Cartesian>& ximesh)
{
    CollisionNodeMulti<Cartesian, volume> node;
    Stencil<Cartesian, volume, double> r1 = createRSwarm<volume>(x1, j1, ximesh);
    Stencil<Cartesian, volume, double> r2 = createRSwarm<volume>(x2, j2, ximesh);
    Stencil<Cartesian, volume, double> w1 = createWSwarm<volume>(x1, j1, ximesh);
    Stencil<Cartesian, volume, double> w2 = createWSwarm<volume>(x2, j2, ximesh);
    node.a  = ximesh.vol(i1) * ximesh.vol(i2) / 
              (interpVol(j1, r1, ximesh) * interpVol(j2, r2, ximesh));
    node.x1 = r1;
    node.x2 = r2;
    node.w1 = w1;
    node.w2 = w2;
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
    Stencil<Cylindrical, volume, double> r1 = createRSwarm<volume>(x1, j1, ximesh);
    Stencil<Cylindrical, volume, double> r2 = createRSwarm<volume>(x2, j2, ximesh);
    node.a  = ximesh.vol(i1) * ximesh.vol(i2) / 
              (interpVol(j1, r1, ximesh) * interpVol(j2, r2, ximesh));
    node.x1 = r1;
    node.x2 = r2;
    return node;
}


template <Volume volume>
struct VolumeConsts {
    static const int rad = 1;
};

template <>
struct VolumeConsts<Grad> {
    static const int rad = 2;
};


template <Volume volume>
std::vector<V3i> makeStencil(V3i vj)
{
    V3i Vi;
    std::vector<V3i> stencil;
    const int rad = VolumeConsts<volume>::rad;
    for (int s1 = -rad; s1 <= rad; ++s1)
        for (int s2 = -rad; s2 <= rad; ++s2)
            for (int s3 = -rad; s3 <= rad; ++s3) {
                const int s = std::abs(s1) + std::abs(s2) + std::abs(s3);
                if ( (s > 0) && (s <= rad) ) 
                    stencil.push_back(vj + V3i(s1, s2, s3));
            }
    return stencil;
}

template <Volume volume>
std::vector<V2i> makeStencil(V2i vj)
{
    V2i Vi;
    std::vector<V2i> stencil;
    const int rad = VolumeConsts<volume>::rad;
    for (int s1 = -rad; s1 <= rad; ++s1)
        for (int s2 = -rad; s2 <= rad; ++s2) {
            const int s = std::abs(s1) + std::abs(s2);
            if ( (s > 0) && (s <= rad) )
                stencil.push_back(vj + V2i(s1, s2));
        }
    return stencil;
}

template <Volume volume, Symmetry symmetry>
std::vector< typename XiMeshMix<symmetry>::Vi >
stencilMulti(typename XiMeshMix<symmetry>::Vi vj, const XiMeshMix<symmetry>& )
{
    typedef typename SymmetryTrait<symmetry>::Vi Vi_simple;
    std::vector<Vi_simple> points = makeStencil<volume>(vj.vi);

    typedef typename XiMeshMix<symmetry>::Vi Vi;
    std::vector<Vi> stencil;

    for (size_t i = 0; i < points.size(); ++i)
        stencil.push_back(Vi(points[i], vj.i));

    return stencil;
}

template <Volume volume, Symmetry symmetry>
std::vector< typename XiMeshRect<symmetry>::Vi >
stencilMulti(typename XiMeshRect<symmetry>::Vi vj, const XiMeshRect<symmetry>& )
{
    return makeStencil<volume>(vj);
}

template <Volume volume, Symmetry symmetry, template <Symmetry> class XiMeshType>
std::tuple< bool, typename CollisionNodeMulti<symmetry, volume>::Si, 
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
        const std::vector<Vi> stencil = stencilMulti<volume, symmetry>(vj, ximesh);
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
        return std::make_tuple(b, j1, x1);
    }
    return std::make_tuple(true, j1, x1);
}

template <Symmetry symmetry> 
double calculateSection(double g, double dot_un, int i1, int i2,
                        const SimpleSection* section,
                        const XiMeshMix<symmetry>& mesh)
{
    double m1 = mesh.m(i1);
    double m2 = mesh.m(i2);
    return g * section->section(g * std::sqrt(2 * m1 * m2 / (m1 + m2)),
                                  dot_un,
                                  mesh.i2ci(i1), mesh.i2ci(i2));
}

template <Symmetry symmetry> 
double calculateSection(double g, double dot_un, int i1, int i2,
                        const SimpleSection* section,
                        const XiMeshRect<symmetry>& mesh)
{
    return g * section->section(g, dot_un);
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
    std::tie(i1, i2, vi, wi, v, w, n) = calcNodeBefore(p, ximesh);

//    std::cout << "v, w = " << v << ' ' << w << std::endl;

    double g;
    V3d    u, o, nn, v2, w2;
    std::tie(u, g, o, nn, v2, w2) = calcNodeCollide(i1, i2, v, w, n, ximesh);

    if (g == 0)
        return SameXi;

//    std::cout << "u, g = " << u << ' ' << g << std::endl;

    Vi vj, wj;
    Vd v1, w1;
    std::tie(vj, wj, v1, w1) = calcNodeAfter(i1, i2, v2, w2, ximesh);

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
    std::tie(b, j1, x1) = searchMulti<volume>(v1, vj, ximesh);
    if (not b) return BadSearch;
    std::tie(b, j2, x2) = searchMulti<volume>(w1, wj, ximesh);
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

    node.c = calculateSection(g, dot(u, n) / g, i1, i2, section, ximesh);
    
    return GoodOne;
}

template <typename CollisionNodeType>
inline double aff(const double ff, const CollisionNodeType n, SymmetryTrait<Cartesian>) {
    return n.a * ff;
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

template <PowMethod powmethod, typename DF, class Nodes, class XiMeshType>
double ciIterMultiCont(DF& f, const Nodes& nodes, const XiMeshType& ximesh)
{
    int i1 = 0, i2 = 0;
    double exclusion = 0;
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
        /*
        const Sd w1 = n.w1;
        const Sd w2 = n.w2;
        */
//        const double e = interpF(n.e1, r1) + interpF(n.e2, r2) - n.ei1 - n.ei2;
//        std::cout << "e = " << e << std::endl;


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
        const Sd y1 = div(f1, n.v1); 
        const Sd y2 = div(f2, n.v2); 

        const Sd ly1 = mylog<powmethod>(y1);
        const Sd ly2 = mylog<powmethod>(y2);

        const double z1 = interpF(ly1, r1);
        const double z2 = interpF(ly2, r2);

        double zz = aff<Node::symmetry>(exp<powmethod>(z1, z2), n);
        zz *=  n.vi1 * n.vi2 / n.a;
        
        double h1 = d * (std::log(fi1/n.vi1) + std::log(fi2/n.vi2) - z1 - z2);

        if (h1 < 0) 
        {
            std::cout << "h1 = " << h1 << std::endl;
            std::cout << "fi1 = " << fi1 << std::endl;
            std::cout << "fi2 = " << fi2 << std::endl;
            std::cout << "g1 =  " << g1 << std::endl;
            std::cout << "g2 =  " << g2 << std::endl;
            std::cout << "r1 =  " << r1 << std::endl;
            std::cout << "r2 =  " << r2 << std::endl;
            std::cout << "n.a = " << n.a << std::endl;
            std::cout << "vi1 = " << n.vi1 << std::endl;
            std::cout << "vi2 = " << n.vi2 << std::endl;
            std::cout << "v1  = " << n.v1 << std::endl;
            std::cout << "v2  = " << n.v2 << std::endl;
            std::cout << "zz, ff1 = " << zz << ' ' << ff1 << std::endl;
        }
*/
/*
        std::cout << "fi1 = " << fi1 << std::endl;
        std::cout << "fi2 = " << fi2 << std::endl;
        std::cout << "d = " << d << std::endl;
*/
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
            /* If this interpolation is bad, then do nothing */
            fIs(f, n.j1, f1);
            fIs(f, n.j2, f2);
            f[n.i1] = fi1;
            f[n.i2] = fi2;

            const double ff1_ = f1.o * f2.o;
            const double d = (ff-ff1_) * n.c;
            exclusion += std::fabs(d);

/*
            const double g1_ = interpF(logf1, r1);
            const double g2_ = interpF(logf2, r2);
            const double ff1_ = aff<Node::symmetry>(exp<powmethod>(g1_, g2_), n);
            d = ( - ff1_ + ff ) * n.c - d;

            if (std::fabs(d) > 1e+8)
            {
                std::cerr << "d: " << d << std::endl;
                std::cerr << "exclusion: " << f[n.i1] << ' ' << f[n.i2] << std::endl;
                std::cerr << "exclusion_proj: " << fAt(f, n.j1) << ' ' << fAt(f, n.j2);
                std::cerr << "r_coeffs: " << r1 << ' ' << r2;
                std::cerr << "w_coeffs: " << n.w1 << ' ' << n.w2;
                std::cerr << "before_proj: " << f1 << ' ' << f2;
                std::cerr << "before: " << fi1 << ' ' << fi2 << std::endl;
                std::cerr << "(ff, ff1_pow/my)*n.c: " << ff*n.c << ' ' << ff1_*n.c << ' ' << ff1*n.c << std::endl;
                std::cerr << "f_alpha', f_beta' (Power): " << std::exp2(g1_) << ' ' << std::exp2(g2_) << std::endl;
                std::cerr << "f_alpha', f_beta'    (My): " << std::exp2(g1) << ' ' << std::exp2(g2) << std::endl;
                std::cerr << "d_pow/my: " << (ff-ff1_)*n.c << ' ' << (ff-ff1)*n.c << std::endl << std::endl;
            }
*/
            ++i2;
        }
        ++i1;
    }
    //std::cout << "i1, i2 = " << i1 << ' ' << i2 << ' ' << exclusion << std::endl;
    //std::cout << "i1, i2 = " << i1 << ' ' << i2 << ' ' << (i2 + 0.0) / i1 << std::endl;
    //return static_cast<double>(i2) / nodes.size();
    return exclusion;
}

#endif
