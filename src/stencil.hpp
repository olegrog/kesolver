#ifndef _STENCIL_HPP_
#define _STENCIL_HPP_

#include <iostream>

#include "symmetry.hpp"
#include "v.hpp"
#include "sse.hpp"

enum PowMethod {
    Std, FastSSE
};

enum Volume {
    Tight,
    Varghese,
    Symmetric,
    Wide,
    Tight2,
    Grad
};

template <Volume volume> struct trait_volume {};

inline std::ostream& operator<<(std::ostream& to, Volume volume) {
    if       (volume == Tight)      return to << "Tight"     ;
    else if  (volume == Symmetric)  return to << "Symmetric" ;
    else if  (volume == Wide)       return to << "Wide"      ;
    else if  (volume == Tight2)     return to << "Tight2"    ;
    else if  (volume == Grad)       return to << "Grad"      ;
    else                            return to << volume      ;
}

template <Symmetry symmetry, Volume volume, typename T>
struct Stencil;

template <typename T>
struct Stencil<Cylindrical, Symmetric, T> {
    T o, x;
    T i, l, j, m;
};

template <typename T>
struct Stencil<Cartesian, Symmetric, T> {
    T o, x;
    T i, l, j, m, k, n;
};

template <typename T>
struct Stencil<Cylindrical, Wide, T> {
    T o, x;
    T i, l, j, m;
    T ij, lj, im, lm;
};

template <typename T>
struct Stencil<Cartesian, Wide, T> {
    T o, x;
    T i, l, j, m, k, n;
    T ij, lj, im, lm;
    T ik, lk, in, ln;
    T jk, mk, jn, mn;
};

template <typename T>
struct Stencil<Cylindrical, Tight, T> {
    T o, x;
    T i, j;
};

template <typename T>
struct Stencil<Cartesian, Tight, T> {
    T o, x;
    T i, j, k, y;
};

template <typename T>
struct Stencil<Cylindrical, Tight2, T> {
    T o, ij;
    T i, l, j, m;
};

template <typename T>
struct Stencil<Cartesian, Tight2, T> {
    T i, l, j, m, k, n;
    T o;
    T ij, ik, jk;
};

template <typename T>
struct Stencil<Cylindrical, Grad, T> {
    T o, ij;
    T i, l, j, m;
    T ii, jj;
};

template <typename T>
struct Stencil<Cartesian, Grad, T> {
    T o, x;
    T i, l, j, m, k, n;
    T ij, ik, jk, ii, jj, kk;
};

template <typename T>
inline std::ostream& operator<<(std::ostream& to, const Stencil<Cartesian, Symmetric, T>& s) {
	to << "[[\n"
       << s.o  << '\n' 
       << s.i  << ' ' << s.j  << ' ' << s.k << ' ' << s.l  << ' ' << s.m  << ' ' << s.n << '\n'
       << "]]\n";
	return to;
}

template <typename T>
inline std::ostream& operator<<(std::ostream& to, const Stencil<Cylindrical, Symmetric, T>& s) {
	to << s.o  << '\n' 
       << s.i  << ' ' << s.j  << ' ' << s.l  << ' ' << s.m  << '\n'
       << std::endl;
	return to;
}

template <class T> inline
std::ostream& operator<<(std::ostream& to, const Stencil<Cartesian, Wide, T>& s) {
	to << "[[\n"
       << s.o  << '\n' 
       << s.i  << ' ' << s.j  << ' ' << s.k << ' ' << s.l  << ' ' << s.m  << ' ' << s.n << '\n'
       << s.ij << ' ' << s.lj << ' ' << s.lm << ' ' << s.im << '\n'
       << s.ik << ' ' << s.lk << ' ' << s.ln << ' ' << s.in << '\n'
       << s.jk << ' ' << s.mk << ' ' << s.mn << ' ' << s.jn << '\n'
       << "]]\n";
	return to;
}

template <class T> inline
std::ostream& operator<<(std::ostream& to, const Stencil<Cylindrical, Wide, T>& s) {
	to << s.o  << '\n' 
       << s.i  << ' ' << s.j  << ' ' << s.l  << ' ' << s.m  << '\n'
       << s.ij << ' ' << s.lj << ' ' << s.lm << ' ' << s.im << '\n'
       << std::endl;
	return to;
}

template <class T> inline
std::ostream& operator<<(std::ostream& to, const Stencil<Cartesian, Tight, T>& s) {
	to << "[[\n"
       << s.o  << s.x << '\n' 
       << s.i  << ' ' << s.j  << ' ' << s.k << '\n'
       << "]]\n";
	return to;
}

template <class T> inline
std::ostream& operator<<(std::ostream& to, const Stencil<Cylindrical, Tight, T>& s) {
	to << s.o  << ' ' << s.x << '\n' 
       << s.i  << ' ' << s.j << '\n'
       << std::endl;
	return to;
}

template <typename T>
inline std::ostream& operator<<(std::ostream& to, const Stencil<Cartesian, Tight2, T>& s) {
	to << "[[\n"
       << s.o  << '\n' 
       << s.i  << ' ' << s.j  << ' ' << s.k << ' ' << s.l  << ' ' << s.m  << ' ' << s.n << '\n'
       << s.ij << ' ' << s.ik << ' ' << s.jk << '\n'
       << "]]\n";
	return to;
}

template <typename T>
inline std::ostream& operator<<(std::ostream& to, const Stencil<Cylindrical, Tight2, T>& s) {
	to << s.o  << '\n' 
       << s.i  << ' ' << s.j  << ' ' << s.l  << ' ' << s.m  << '\n'
       << s.ij << '\n'
       << std::endl;
	return to;
}

template <typename T>
inline std::ostream& operator<<(std::ostream& to, const Stencil<Cartesian, Grad, T>& s) {
	to << "[[\n"
       << s.o  << '\n' 
       << s.i  << ' ' << s.j  << ' ' << s.k << ' ' << s.l  << ' ' << s.m  << ' ' << s.n << '\n'
       << s.ij << ' ' << s.ik << ' ' << s.jk << '\n'
       << s.ii << ' ' << s.jj << ' ' << s.kk << '\n'
       << "]]\n";
	return to;
}

template <typename T>
inline std::ostream& operator<<(std::ostream& to, const Stencil<Cylindrical, Grad, T>& s) {
	to << s.o  << '\n' 
       << s.i  << ' ' << s.j  << ' ' << s.l  << ' ' << s.m  << '\n'
       << s.ij << '\n'
       << s.ii << ' ' << s.jj << '\n'
       << std::endl;
	return to;
}

template <Volume volume, Symmetry symmetry, template <Symmetry> class XiMeshType>
Stencil<symmetry, volume, int>
makeISwarm(const typename XiMeshType<symmetry>::Vi xi,
           const typename XiMeshType<symmetry>::Vd x,
           const XiMeshType<symmetry>& ximesh)
{
    return makeISwarm(xi, x, ximesh, trait_volume<volume>());
}

template <template <Symmetry symmetry> class XiMeshType>
Stencil<Cylindrical, Symmetric, int>
makeISwarm(const typename XiMeshType<Cylindrical>::Vi xi,
           const V2d x,
           const XiMeshType<Cylindrical>& ximesh,
           trait_volume<Symmetric> )
{
    Stencil<Cylindrical, Symmetric, int> i;

    i.o  = ximesh(xi);

    i.i  = ximesh(xi + V2i( 1,  0));
    i.l  = ximesh(xi - V2i( 1,  0));
    i.j  = ximesh(xi + V2i( 0,  1));
    i.m  = ximesh(xi - V2i( 0,  1));
    
//    std::cout << "i = " << i << std::endl;

    return i;
}

template <template <Symmetry symmetry> class XiMeshType>
Stencil<Cartesian, Symmetric, int>
makeISwarm(const typename XiMeshType<Cartesian>::Vi xi,
           const V3d x,
           const XiMeshType<Cartesian>& ximesh,
           trait_volume<Symmetric> )
{
    Stencil<Cartesian, Symmetric, int> i;

    i.o  = ximesh(xi);

    i.i  = ximesh(xi + V3i( 1,  0,  0));
    i.l  = ximesh(xi - V3i( 1,  0,  0));
    i.j  = ximesh(xi + V3i( 0,  1,  0));
    i.m  = ximesh(xi - V3i( 0,  1,  0));
    i.k  = ximesh(xi + V3i( 0,  0,  1));
    i.n  = ximesh(xi - V3i( 0,  0,  1));

    return i;
}

template <template <Symmetry symmetry> class XiMeshType>
Stencil<Cylindrical, Wide, int>
makeISwarm(const typename XiMeshType<Cylindrical>::Vi xi,
           const V2d x,
           const XiMeshType<Cylindrical>& ximesh,
           trait_volume<Wide> )
{
    Stencil<Cylindrical, Wide, int> i;

    i.o  = ximesh(xi);

    i.i  = ximesh(xi + V2i( 1,  0));
    i.l  = ximesh(xi - V2i( 1,  0));
    i.j  = ximesh(xi + V2i( 0,  1));
    i.m  = ximesh(xi - V2i( 0,  1));

    i.ij = ximesh(xi + V2i( 1,  1));
    i.lj = ximesh(xi + V2i(-1,  1));
    i.lm = ximesh(xi + V2i(-1, -1));
    i.im = ximesh(xi + V2i( 1, -1));

    return i;
}

template <template <Symmetry symmetry> class XiMeshType>
Stencil<Cartesian, Wide, int>
makeISwarm(const typename XiMeshType<Cartesian>::Vi xi,
           const V3d x,
           const XiMeshType<Cartesian>& ximesh,
           trait_volume<Wide> )
{
    Stencil<Cartesian, Wide, int> i;

    i.o  = ximesh(xi);

    i.i  = ximesh(xi + V3i( 1,  0,  0));
    i.l  = ximesh(xi - V3i( 1,  0,  0));
    i.j  = ximesh(xi + V3i( 0,  1,  0));
    i.m  = ximesh(xi - V3i( 0,  1,  0));
    i.k  = ximesh(xi + V3i( 0,  0,  1));
    i.n  = ximesh(xi - V3i( 0,  0,  1));

    i.ij = ximesh(xi + V3i( 1,  1,  0));
    i.lj = ximesh(xi + V3i(-1,  1,  0));
    i.lm = ximesh(xi + V3i(-1, -1,  0));
    i.im = ximesh(xi + V3i( 1, -1,  0));

    i.ik = ximesh(xi + V3i( 1,  0,  1));
    i.lk = ximesh(xi + V3i(-1,  0,  1));
    i.ln = ximesh(xi + V3i(-1,  0, -1));
    i.in = ximesh(xi + V3i( 1,  0, -1));

    i.jk = ximesh(xi + V3i( 0,  1,  1));
    i.mk = ximesh(xi + V3i( 0, -1,  1));
    i.mn = ximesh(xi + V3i( 0, -1, -1));
    i.jn = ximesh(xi + V3i( 0,  1, -1));

    return i;
}

template <template <Symmetry symmetry> class XiMeshType>
Stencil<Cylindrical, Tight, int>
makeISwarm(const typename XiMeshType<Cylindrical>::Vi xi,
           const V2d x,
           const XiMeshType<Cylindrical>& ximesh,
           trait_volume<Tight> )
{
    Stencil<Cylindrical, Tight, int> i;
    V2i d = sign(x);

    i.o  = ximesh(xi);
    i.x  = ximesh(xi - d);

    i.i  = ximesh(xi + V2i( d[0],  0    ));
    i.j  = ximesh(xi + V2i( 0,     d[1] ));

    return i;
}

template <template <Symmetry symmetry> class XiMeshType>
Stencil<Cartesian, Tight, int>
makeISwarm(const typename XiMeshType<Cartesian>::Vi xi,
           const V3d x,
           const XiMeshType<Cartesian>& ximesh,
           trait_volume<Tight> )
{
    Stencil<Cartesian, Tight, int> i;
    V3i d = sign(x);

    i.o  = ximesh(xi);
    i.x  = ximesh(xi - d);

    i.i  = ximesh(xi + V3i( d[0],  0,  0));
    i.j  = ximesh(xi + V3i( 0,  d[1],  0));
    i.k  = ximesh(xi + V3i( 0,  0,  d[2]));

    return i;
}

template <template <Symmetry symmetry> class XiMeshType>
Stencil<Cylindrical, Tight2, int>
makeISwarm(const typename XiMeshType<Cylindrical>::Vi xi,
           const V2d x,
           const XiMeshType<Cylindrical>& ximesh,
           trait_volume<Tight2> )
{
    Stencil<Cylindrical, Tight2, int> i;
    V2i d = sign(x);

    i.o  = ximesh(xi);
    i.ij = ximesh(xi + d);

    i.i  = ximesh(xi + V2i( d[0],  0    ));
    i.j  = ximesh(xi + V2i( 0,     d[1] ));
    i.l  = ximesh(xi + V2i( - d[0],  0    ));
    i.m  = ximesh(xi + V2i( 0,     - d[1] ));

    return i;
}

template <template <Symmetry symmetry> class XiMeshType>
Stencil<Cartesian, Tight2, int>
makeISwarm(const typename XiMeshType<Cartesian>::Vi xi,
           const V3d x,
           const XiMeshType<Cartesian>& ximesh,
           trait_volume<Tight2> )
{
    Stencil<Cartesian, Tight2, int> i;
    V3i d = sign(x);

    i.o  = ximesh(xi);

    i.i  = ximesh(xi + V3i( d[0],  0,  0));
    i.j  = ximesh(xi + V3i( 0,  d[1],  0));
    i.k  = ximesh(xi + V3i( 0,  0,  d[2]));

    i.l  = ximesh(xi + V3i( - d[0],  0,  0));
    i.m  = ximesh(xi + V3i( 0,  - d[1],  0));
    i.n  = ximesh(xi + V3i( 0,  0,  - d[2]));

    i.ij = ximesh(xi + V3i( d[0],  d[1],  0));
    i.ik = ximesh(xi + V3i( d[0],  0,  d[2]));
    i.jk = ximesh(xi + V3i( 0,  d[1],  d[2]));

    return i;
}

template <template <Symmetry symmetry> class XiMeshType>
Stencil<Cylindrical, Grad, int>
makeISwarm(const typename XiMeshType<Cylindrical>::Vi xi,
           const V2d x,
           const XiMeshType<Cylindrical>& ximesh,
           trait_volume<Grad> )
{
    Stencil<Cylindrical, Grad, int> i;
    V2i d = sign(x);

    i.o  = ximesh(xi);
    i.ij = ximesh(xi + d);

    i.i  = ximesh(xi + V2i( d[0],  0    ));
    i.j  = ximesh(xi + V2i( 0,     d[1] ));
    i.l  = ximesh(xi + V2i( - d[0],  0    ));
    i.m  = ximesh(xi + V2i( 0,     - d[1] ));

    i.ii = ximesh(xi + V2i( 2*d[0],  0      ));
    i.jj = ximesh(xi + V2i( 0,       2*d[1] ));

    return i;
}

template <template <Symmetry symmetry> class XiMeshType>
Stencil<Cartesian, Grad, int>
makeISwarm(const typename XiMeshType<Cartesian>::Vi xi,
           const V3d x,
           const XiMeshType<Cartesian>& ximesh,
           trait_volume<Grad> )
{
    Stencil<Cartesian, Grad, int> i;
    V3i d = sign(x);

    i.o  = ximesh(xi);

    i.i  = ximesh(xi + V3i( d[0],  0,  0));
    i.j  = ximesh(xi + V3i( 0,  d[1],  0));
    i.k  = ximesh(xi + V3i( 0,  0,  d[2]));

    i.l  = ximesh(xi + V3i( - d[0],  0,  0));
    i.m  = ximesh(xi + V3i( 0,  - d[1],  0));
    i.n  = ximesh(xi + V3i( 0,  0,  - d[2]));

    i.ij = ximesh(xi + V3i( d[0],  d[1],  0));
    i.ik = ximesh(xi + V3i( d[0],  0,  d[2]));
    i.jk = ximesh(xi + V3i( 0,  d[1],  d[2]));

    i.ii = ximesh(xi + V3i( 2*d[0],  0,  0));
    i.jj = ximesh(xi + V3i( 0,  2*d[1],  0));
    i.kk = ximesh(xi + V3i( 0,  0,  2*d[2]));

    return i;
}

template <Volume volume>
Stencil<Cylindrical, volume, double> makeRSwarm(const V2d x);

template <Volume volume>
Stencil<Cartesian, volume, double> makeRSwarm(const V3d x);

template <>
Stencil<Cylindrical, Symmetric, double> makeRSwarm<Symmetric>(const V2d x) {
    Stencil<Cylindrical, Symmetric, double> s;
    s.o  =  1 - sqr(x);

    s.i  =  0.5 * x[0] * (x[0] + 1);
    s.l  =  0.5 * x[0] * (x[0] - 1);
    s.j  =  0.5 * x[1] * (x[1] + 1);
    s.m  =  0.5 * x[1] * (x[1] - 1);

    return s;
}

template <>
Stencil<Cartesian, Symmetric, double> makeRSwarm<Symmetric>(const V3d x) {
    Stencil<Cartesian, Symmetric, double> s;
    s.o  =  1 - sqr(x);

    s.i  =  0.5 * x[0] * (x[0] + 1);
    s.l  =  0.5 * x[0] * (x[0] - 1);
    s.j  =  0.5 * x[1] * (x[1] + 1);
    s.m  =  0.5 * x[1] * (x[1] - 1);
    s.k  =  0.5 * x[2] * (x[2] + 1);
    s.n  =  0.5 * x[2] * (x[2] - 1);

    return s;
}

template <>
Stencil<Cylindrical, Wide, double> makeRSwarm<Wide>(const V2d x) {
    Stencil<Cylindrical, Wide, double> s;
    s.o  =  1 - sqr(x);

    s.i  =  0.5 * x[0] * (x[0] + 1);
    s.l  =  0.5 * x[0] * (x[0] - 1);
    s.j  =  0.5 * x[1] * (x[1] + 1);
    s.m  =  0.5 * x[1] * (x[1] - 1);

    s.ij =  0.25 * x[0] * x[1];
    s.lj = -0.25 * x[0] * x[1];
    s.lm =  0.25 * x[0] * x[1];
    s.im = -0.25 * x[0] * x[1];

    return s;
}

template <>
Stencil<Cartesian, Wide, double> makeRSwarm<Wide>(const V3d x) {
    Stencil<Cartesian, Wide, double> s;
    s.o  =  1 - sqr(x);

    s.i  =  0.5 * x[0] * (x[0] + 1);
    s.l  =  0.5 * x[0] * (x[0] - 1);
    s.j  =  0.5 * x[1] * (x[1] + 1);
    s.m  =  0.5 * x[1] * (x[1] - 1);
    s.k  =  0.5 * x[2] * (x[2] + 1);
    s.n  =  0.5 * x[2] * (x[2] - 1);

    s.ij =  0.25 * x[0] * x[1];
    s.lj = -0.25 * x[0] * x[1];
    s.lm =  0.25 * x[0] * x[1];
    s.im = -0.25 * x[0] * x[1];

    s.ik =  0.25 * x[0] * x[2];
    s.lk = -0.25 * x[0] * x[2];
    s.ln =  0.25 * x[0] * x[2];
    s.in = -0.25 * x[0] * x[2];

    s.jk =  0.25 * x[2] * x[1];
    s.mk = -0.25 * x[2] * x[1];
    s.mn =  0.25 * x[2] * x[1];
    s.jn = -0.25 * x[2] * x[1];

    return s;
}

template <>
Stencil<Cylindrical, Tight, double> makeRSwarm<Tight>(const V2d x) {
    Stencil<Cylindrical, Tight, double> s;
    const double square = sqr(x);
    const double sum    = norm_abs(x);
    const double rext   = - 1. / 4. * (sum - square);

    s.o  =  1 - sum - 3*rext;
    s.x  =  rext;

    s.i  =  templ_abs(x[0]) + rext;
    s.j  =  templ_abs(x[1]) + rext;

//    std::cout << "x, s = " << x << ' ' << s << std::endl;

    return s;
}

template <>
Stencil<Cartesian, Tight, double> makeRSwarm<Tight>(const V3d x) {
    Stencil<Cartesian, Tight, double> s;
    const double square = sqr(x);
    const double sum    = norm_abs(x);
    const double rext   = - 1. / 6. * (sum - square);

    s.o  =  1 - sum - 4*rext;
    s.x  =  rext;

    s.i  =  templ_abs(x[0]) + rext;
    s.j  =  templ_abs(x[1]) + rext;
    s.k  =  templ_abs(x[2]) + rext;

    return s;
}

template <>
Stencil<Cylindrical, Tight2, double> makeRSwarm<Tight2>(const V2d r) {
    Stencil<Cylindrical, Tight2, double> s;
    
    double x = templ_abs(r[0]);
    double y = templ_abs(r[1]);

    double xy = x * y;

    s.o  =  1 - sqr(r) + xy;

    s.i  =  0.5 * x * (x + 1) - xy;
    s.l  =  0.5 * x * (x - 1);
    s.j  =  0.5 * y * (y + 1) - xy;
    s.m  =  0.5 * y * (y - 1);

    s.ij =  xy;

    return s;
}

template <>
Stencil<Cartesian, Tight2, double> makeRSwarm<Tight2>(const V3d r) {
    Stencil<Cartesian, Tight2, double> s;
    
    double x = templ_abs(r[0]);
    double y = templ_abs(r[1]);
    double z = templ_abs(r[2]);

    double xy = x * y;
    double xz = x * z;
    double yz = y * z;

    s.o  =  1 - sqr(r) + xy + xz + yz;

    s.i  =  0.5 * x * (x + 1) - (xy + xz);
    s.l  =  0.5 * x * (x - 1);
    s.j  =  0.5 * y * (y + 1) - (xy + yz);
    s.m  =  0.5 * y * (y - 1);
    s.k  =  0.5 * z * (z + 1) - (xz + yz);
    s.n  =  0.5 * z * (z - 1);

    s.ij =  xy;
    s.ik =  xz;
    s.jk =  yz;

    return s;
}

template <>
Stencil<Cylindrical, Grad, double> makeRSwarm<Grad>(const V2d r) {
    Stencil<Cylindrical, Grad, double> s;
    
    double x = templ_abs(r[0]);
    double y = templ_abs(r[1]);

    double xy = x * y;

    double qx = x * ( 1 - x*x );
    double qy = y * ( 1 - y*y );

    s.o  =  1 - sqr(r) + xy - 0.5 * (qx + qy);

    s.i  =  0.5 * x * (x + 1) - xy + qx / 2;
    s.l  =  0.5 * x * (x - 1)      + qx / 6;
    s.j  =  0.5 * y * (y + 1) - xy + qy / 2;
    s.m  =  0.5 * y * (y - 1)      + qy / 6;

    s.ij =  xy;

    s.ii =  - qx / 6.;
    s.jj =  - qy / 6.;

    return s;
}

template <>
Stencil<Cartesian, Grad, double> makeRSwarm<Grad>(const V3d r) {
    Stencil<Cartesian, Grad, double> s;
    
    double x = templ_abs(r[0]);
    double y = templ_abs(r[1]);
    double z = templ_abs(r[2]);

    double xy = x * y;
    double xz = x * z;
    double yz = y * z;

    double qx = x * ( 1 - x*x );
    double qy = y * ( 1 - y*y );
    double qz = z * ( 1 - z*z );

    s.o  =  1 - sqr(r) + xy + xz + yz - 0.5 * (qx + qy + qz);

    s.i  =  0.5 * x * (x + 1) - (xy + xz) + qx / 2;
    s.l  =  0.5 * x * (x - 1)             + qx / 6;
    s.j  =  0.5 * y * (y + 1) - (xy + yz) + qy / 2;
    s.m  =  0.5 * y * (y - 1)             + qy / 6;
    s.k  =  0.5 * z * (z + 1) - (xz + yz) + qz / 2;
    s.n  =  0.5 * z * (z - 1)             + qz / 6;

    s.ij =  xy;
    s.ik =  xz;
    s.jk =  yz;

    s.ii =  - qx / 6.;
    s.jj =  - qy / 6.;
    s.kk =  - qz / 6.;

    return s;
}

template <typename T>
inline bool lessZero(const Stencil<Cylindrical, Symmetric, T> i) {
    return ((i.o < 0)  ||
            (i.i < 0)  || (i.j < 0)  || 
            (i.l < 0)  || (i.m < 0)  );
}

template <typename T>
inline bool lessZero(const Stencil<Cartesian, Symmetric, T> i) {
    return ((i.o < 0)  ||
            (i.i < 0)  || (i.j < 0)  || (i.k < 0)  || 
            (i.l < 0)  || (i.m < 0)  || (i.n < 0)  );
}

template <typename T>
inline bool lessZero(const Stencil<Cylindrical, Wide, T> i) {
    return ((i.o < 0)  ||
            (i.i < 0)  || (i.j < 0)  || 
            (i.l < 0)  || (i.m < 0)  ||
            (i.ij < 0) || (i.lj < 0) || (i.lm < 0) || (i.im < 0));
}

template <typename T>
inline bool lessZero(const Stencil<Cartesian, Wide, T> i) {
    return ((i.o < 0)  ||
            (i.i < 0)  || (i.j < 0)  || (i.k < 0)  || 
            (i.l < 0)  || (i.m < 0)  || (i.n < 0)  ||
            (i.ij < 0) || (i.lj < 0) || (i.lm < 0) || (i.im < 0) ||
            (i.ik < 0) || (i.lk < 0) || (i.ln < 0) || (i.in < 0) ||
            (i.jk < 0) || (i.mk < 0) || (i.mn < 0) || (i.jn < 0));
}

template <typename T>
inline bool lessZero(const Stencil<Cylindrical, Tight, T> i) {
    return ((i.o < 0)  || (i.x < 0)  ||
            (i.i < 0)  || (i.j < 0)  );
}

template <typename T>
inline bool lessZero(const Stencil<Cartesian, Tight, T> i) {
    return ((i.o < 0)  || (i.x < 0)  ||
            (i.i < 0)  || (i.j < 0)  || (i.k < 0) );
}

template <typename T>
inline bool lessZero(const Stencil<Cylindrical, Tight2, T> i) {
    return ((i.o < 0)  ||
            (i.i < 0)  || (i.j < 0)  || 
            (i.l < 0)  || (i.m < 0)  ||
            (i.ij < 0) );
}

template <typename T>
inline bool lessZero(const Stencil<Cartesian, Tight2, T> i) {
    return ((i.o < 0)  ||
            (i.i < 0)  || (i.j < 0)  || (i.k < 0)  || 
            (i.l < 0)  || (i.m < 0)  || (i.n < 0)  ||
            (i.ij < 0) || (i.ik < 0) || (i.jk < 0) );
}

template <typename T>
inline bool lessZero(const Stencil<Cylindrical, Grad, T> i) {
    return ((i.o < 0)  ||
            (i.i < 0)  || (i.j < 0)  || 
            (i.l < 0)  || (i.m < 0)  ||
            (i.ii < 0) || (i.jj < 0) || 
            (i.ij < 0) );
}

template <typename T>
inline bool lessZero(const Stencil<Cartesian, Grad, T> i) {
    return ((i.o < 0)  ||
            (i.i < 0)  || (i.j < 0)  || (i.k < 0)  || 
            (i.l < 0)  || (i.m < 0)  || (i.n < 0)  ||
            (i.ii < 0) || (i.jj < 0) || (i.kk < 0) || 
            (i.ij < 0) || (i.ik < 0) || (i.jk < 0) );
}

template <typename F>
inline Stencil<Cylindrical, Symmetric, double>
fAt(const F& f, const Stencil<Cylindrical, Symmetric, int> j)
{
    Stencil<Cylindrical, Symmetric, double> s;

    s.o = f[j.o];
    s.i = f[j.i];
    s.j = f[j.j];
    s.l = f[j.l];
    s.m = f[j.m];

    return s;
}

template <typename F>
inline Stencil<Cartesian, Symmetric, double>
fAt(const F& f, const Stencil<Cartesian, Symmetric, int> j)
{
    Stencil<Cartesian, Symmetric, double> s;

    s.o = f[j.o];
    s.i = f[j.i];
    s.j = f[j.j];
    s.k = f[j.k];
    s.l = f[j.l];
    s.m = f[j.m];
    s.n = f[j.n];

    return s;
}

template <typename F>
inline Stencil<Cylindrical, Wide, double>
fAt(const F& f, const Stencil<Cylindrical, Wide, int> j)
{
    Stencil<Cylindrical, Wide, double> s;

    s.o = f[j.o];
    s.i = f[j.i];
    s.j = f[j.j];
    s.l = f[j.l];
    s.m = f[j.m];

    s.ij = f[j.ij];
    s.lj = f[j.lj];
    s.lm = f[j.lm];
    s.im = f[j.im];

    return s;
}

template <typename F>
inline Stencil<Cartesian, Wide, double>
fAt(const F& f, const Stencil<Cartesian, Wide, int> j)
{
    Stencil<Cartesian, Wide, double> s;

    s.o = f[j.o];
    s.i = f[j.i];
    s.j = f[j.j];
    s.k = f[j.k];
    s.l = f[j.l];
    s.m = f[j.m];
    s.n = f[j.n];

    s.ij = f[j.ij];
    s.lj = f[j.lj];
    s.lm = f[j.lm];
    s.im = f[j.im];
                 
    s.ik = f[j.ik];
    s.lk = f[j.lk];
    s.ln = f[j.ln];
    s.in = f[j.in];
                 
    s.jk = f[j.jk];
    s.mk = f[j.mk];
    s.mn = f[j.mn];
    s.jn = f[j.jn];

    return s;
}

template <typename F>
inline Stencil<Cylindrical, Tight, double>
fAt(const F& f, const Stencil<Cylindrical, Tight, int> j)
{
    Stencil<Cylindrical, Tight, double> s;

    s.o = f[j.o];
    s.x = f[j.x];
    s.i = f[j.i];
    s.j = f[j.j];

    return s;
}

template <typename F>
inline Stencil<Cartesian, Tight, double>
fAt(const F& f, const Stencil<Cartesian, Tight, int> j)
{
    Stencil<Cartesian, Tight, double> s;

    s.o = f[j.o];
    s.x = f[j.x];
    s.i = f[j.i];
    s.j = f[j.j];
    s.k = f[j.k];

    return s;
}

template <typename F>
inline Stencil<Cylindrical, Tight2, double>
fAt(const F& f, const Stencil<Cylindrical, Tight2, int> j)
{
    Stencil<Cylindrical, Tight2, double> s;

    s.o = f[j.o];
    s.ij = f[j.ij];
    s.i = f[j.i];
    s.j = f[j.j];
    s.l = f[j.l];
    s.m = f[j.m];

    return s;
}

template <typename F>
inline Stencil<Cartesian, Tight2, double>
fAt(const F& f, const Stencil<Cartesian, Tight2, int> j)
{
    Stencil<Cartesian, Tight2, double> s;

    s.i = f[j.i];
    s.j = f[j.j];
    s.k = f[j.k];
    s.l = f[j.l];
    s.m = f[j.m];
    s.n = f[j.n];
    s.o = f[j.o];
    s.ij = f[j.ij];
    s.ik = f[j.ik];
    s.jk = f[j.jk];

    return s;
}

template <typename F>
inline Stencil<Cylindrical, Grad, double>
fAt(const F& f, const Stencil<Cylindrical, Grad, int> j)
{
    Stencil<Cylindrical, Grad, double> s;

    s.o = f[j.o];
    s.ij = f[j.ij];
    s.i = f[j.i];
    s.j = f[j.j];
    s.l = f[j.l];
    s.m = f[j.m];
    s.ii = f[j.ii];
    s.jj = f[j.jj];

    return s;
}

template <typename F>
inline Stencil<Cartesian, Grad, double>
fAt(const F& f, const Stencil<Cartesian, Grad, int> j)
{
    Stencil<Cartesian, Grad, double> s;

    s.i = f[j.i];
    s.j = f[j.j];
    s.k = f[j.k];
    s.l = f[j.l];
    s.m = f[j.m];
    s.n = f[j.n];
    s.o = f[j.o];
    s.ij = f[j.ij];
    s.ik = f[j.ik];
    s.jk = f[j.jk];
    s.ii = f[j.ii];
    s.jj = f[j.jj];
    s.kk = f[j.kk];

    return s;
}


template <typename F, typename T>
inline void fIs(F& f, const Stencil<Cylindrical, Symmetric, int> j,
                      const Stencil<Cylindrical, Symmetric, T>   s)
{
    f[j.o] = s.o;
    f[j.i] = s.i;
    f[j.j] = s.j;
    f[j.l] = s.l;
    f[j.m] = s.m;
}

template <typename F, typename T>
inline void fIs(F& f, const Stencil<Cartesian, Symmetric, int> j,
                      const Stencil<Cartesian, Symmetric, T>   s)
{
    f[j.o] = s.o;
    f[j.i] = s.i;
    f[j.j] = s.j;
    f[j.k] = s.k;
    f[j.l] = s.l;
    f[j.m] = s.m;
    f[j.n] = s.n;
}

template <typename F, typename T>
inline void fIs(F& f, const Stencil<Cylindrical, Wide, int> j,
                      const Stencil<Cylindrical, Wide, T>   s)
{
    f[j.o] = s.o;
    f[j.i] = s.i;
    f[j.j] = s.j;
    f[j.l] = s.l;
    f[j.m] = s.m;

    f[j.ij] = s.ij;
    f[j.lj] = s.lj;
    f[j.lm] = s.lm;
    f[j.im] = s.im;
}

template <typename F, typename T>
inline void fIs(F& f, const Stencil<Cartesian, Wide, int> j,
                      const Stencil<Cartesian, Wide, T>   s)
{
    f[j.o] = s.o;
    f[j.i] = s.i;
    f[j.j] = s.j;
    f[j.k] = s.k;
    f[j.l] = s.l;
    f[j.m] = s.m;
    f[j.n] = s.n;

    f[j.ij] = s.ij;
    f[j.lj] = s.lj;
    f[j.lm] = s.lm;
    f[j.im] = s.im;
                  
    f[j.ik] = s.ik;
    f[j.lk] = s.lk;
    f[j.ln] = s.ln;
    f[j.in] = s.in;
                  
    f[j.jk] = s.jk;
    f[j.mk] = s.mk;
    f[j.mn] = s.mn;
    f[j.jn] = s.jn;
}

template <typename F, typename T>
inline void fIs(F& f, const Stencil<Cylindrical, Tight, int> j,
                      const Stencil<Cylindrical, Tight, T>   s)
{
    f[j.o] = s.o;
    f[j.x] = s.x;
    f[j.i] = s.i;
    f[j.j] = s.j;
}

template <typename F, typename T>
inline void fIs(F& f, const Stencil<Cartesian, Tight, int> j,
                      const Stencil<Cartesian, Tight, T>   s)
{
    f[j.o] = s.o;
    f[j.x] = s.x;
    f[j.i] = s.i;
    f[j.j] = s.j;
    f[j.k] = s.k;
}

template <typename F, typename T>
inline void fIs(F& f, const Stencil<Cylindrical, Tight2, int> j,
                      const Stencil<Cylindrical, Tight2, T>   s)
{
    f[j.o]  = s.o;
    f[j.ij] = s.ij;
    f[j.i]  = s.i;
    f[j.j]  = s.j;
    f[j.l]  = s.l;
    f[j.m]  = s.m;
}

template <typename F, typename T>
inline void fIs(F& f, const Stencil<Cartesian, Tight2, int> j,
                      const Stencil<Cartesian, Tight2, T>   s)
{
    f[j.i] = s.i;
    f[j.j] = s.j;
    f[j.k] = s.k;
    f[j.l] = s.l;
    f[j.m] = s.m;
    f[j.n] = s.n;
    f[j.o] = s.o;
    f[j.ij] = s.ij;
    f[j.ik] = s.ik;
    f[j.jk] = s.jk;
}

template <typename F, typename T>
inline void fIs(F& f, const Stencil<Cylindrical, Grad, int> j,
                      const Stencil<Cylindrical, Grad, T>   s)
{
    f[j.i]  = s.i;
    f[j.j]  = s.j;
    f[j.l]  = s.l;
    f[j.m]  = s.m;
    f[j.o]  = s.o;
    f[j.ij] = s.ij;
    f[j.ii] = s.ii;
    f[j.jj] = s.jj;
}

template <typename F, typename T>
inline void fIs(F& f, const Stencil<Cartesian, Grad, int> j,
                      const Stencil<Cartesian, Grad, T>   s)
{
    f[j.o] = s.o;
    f[j.i] = s.i;
    f[j.j] = s.j;
    f[j.k] = s.k;
    f[j.l] = s.l;
    f[j.m] = s.m;
    f[j.n] = s.n;
    f[j.ij] = s.ij;
    f[j.ik] = s.ik;
    f[j.jk] = s.jk;
    f[j.ii] = s.ii;
    f[j.jj] = s.jj;
    f[j.kk] = s.kk;
}

inline double interpF(const Stencil<Cylindrical, Symmetric, double> f,
                      const Stencil<Cylindrical, Symmetric, double> r) {
    return  r.o  * f.o  + 
            r.i  * f.i  + 
            r.j  * f.j  + 
            r.l  * f.l  + 
            r.m  * f.m  ;
}

inline double interpF(const Stencil<Cartesian, Symmetric, double> f,
                      const Stencil<Cartesian, Symmetric, double> r) {
    return  r.o  * f.o  + 
            r.i  * f.i  + 
            r.j  * f.j  + 
            r.k  * f.k  + 
            r.l  * f.l  + 
            r.m  * f.m  + 
            r.n  * f.n  ; 
}

inline double interpF(const Stencil<Cylindrical, Wide, double> f,
                      const Stencil<Cylindrical, Wide, double> r) {
    return  r.o  * f.o  + 
            r.i  * f.i  + 
            r.j  * f.j  + 
            r.l  * f.l  + 
            r.m  * f.m  + 

            r.ij * f.ij + 
            r.lj * f.lj + 
            r.lm * f.lm + 
            r.im * f.im ;
}

inline double interpF(const Stencil<Cartesian, Wide, double> f,
                      const Stencil<Cartesian, Wide, double> r) {
    return  r.o  * f.o  + 
            r.i  * f.i  + 
            r.j  * f.j  + 
            r.k  * f.k  + 
            r.l  * f.l  + 
            r.m  * f.m  + 
            r.n  * f.n  + 

            r.ij * f.ij + 
            r.lj * f.lj + 
            r.lm * f.lm + 
            r.im * f.im + 

            r.ik * f.ik + 
            r.lk * f.lk + 
            r.ln * f.ln + 
            r.in * f.in + 

            r.jk * f.jk + 
            r.mk * f.mk + 
            r.mn * f.mn + 
            r.jn * f.jn ;
}

inline double interpF(const Stencil<Cylindrical, Tight, double> f,
                      const Stencil<Cylindrical, Tight, double> r) {
    return  r.o  * f.o  + 
            r.x  * f.x  + 
            r.i  * f.i  + 
            r.j  * f.j  ;
}

inline double interpF(const Stencil<Cartesian, Tight, double> f,
                      const Stencil<Cartesian, Tight, double> r) {
    return  r.o  * f.o  + 
            r.x  * f.x  + 
            r.i  * f.i  + 
            r.j  * f.j  + 
            r.k  * f.k  ;
}

inline double interpF(const Stencil<Cylindrical, Tight2, double> f,
                      const Stencil<Cylindrical, Tight2, double> r) {
    return  r.o  * f.o  + 
            r.i  * f.i  + 
            r.j  * f.j  + 
            r.l  * f.l  + 
            r.m  * f.m  +
            r.ij * f.ij ;
}

inline double interpF(const Stencil<Cartesian, Tight2, double> f,
                      const Stencil<Cartesian, Tight2, double> r) {
    return  r.o  * f.o  + 
            r.i  * f.i  + 
            r.j  * f.j  + 
            r.k  * f.k  + 
            r.l  * f.l  + 
            r.m  * f.m  + 
            r.n  * f.n  + 
            r.ij * f.ij + 
            r.ik * f.ik + 
            r.jk * f.jk ; 
}

inline double interpF(const Stencil<Cylindrical, Grad, double> f,
                      const Stencil<Cylindrical, Grad, double> r) {
    return  r.o  * f.o  + 
            r.i  * f.i  + 
            r.j  * f.j  + 
            r.l  * f.l  + 
            r.m  * f.m  +
            r.ij * f.ij +
            r.ii * f.ii +
            r.jj * f.jj ;
}

inline double interpF(const Stencil<Cartesian, Grad, double> f,
                      const Stencil<Cartesian, Grad, double> r) {
    return  r.o  * f.o  + 
            r.i  * f.i  + 
            r.j  * f.j  + 
            r.k  * f.k  + 
            r.l  * f.l  + 
            r.m  * f.m  + 
            r.n  * f.n  + 
            r.ij * f.ij + 
            r.ik * f.ik + 
            r.jk * f.jk +
            r.ii * f.ii +
            r.jj * f.jj +
            r.kk * f.kk ; 
}

template <typename XiMesh>
inline double interpVol(const Stencil<Cylindrical, Symmetric, int>    j,
                        const Stencil<Cylindrical, Symmetric, double> r,
                        const XiMesh& ximesh) {
    return  std::exp(
                r.o  * ximesh.volLog(j.o ) + 

                r.i  * ximesh.volLog(j.i ) + 
                r.j  * ximesh.volLog(j.j ) + 
                r.l  * ximesh.volLog(j.l ) + 
                r.m  * ximesh.volLog(j.m )
            );
}

template <typename XiMesh>
inline double interpVol(const Stencil<Cartesian, Symmetric, int>    j,
                        const Stencil<Cartesian, Symmetric, double> r,
                        const XiMesh& ximesh) {
    return  std::exp(
                r.o  * ximesh.volLog(j.o ) + 

                r.i  * ximesh.volLog(j.i ) + 
                r.j  * ximesh.volLog(j.j ) + 
                r.k  * ximesh.volLog(j.k ) + 
                r.l  * ximesh.volLog(j.l ) + 
                r.m  * ximesh.volLog(j.m ) + 
                r.n  * ximesh.volLog(j.n )  
            );
}

template <typename XiMesh>
inline double interpVol(const Stencil<Cylindrical, Wide, int>    j,
                        const Stencil<Cylindrical, Wide, double> r,
                        const XiMesh& ximesh) {
    return  std::exp(
                r.o  * ximesh.volLog(j.o ) + 

                r.i  * ximesh.volLog(j.i ) + 
                r.j  * ximesh.volLog(j.j ) + 
                r.l  * ximesh.volLog(j.l ) + 
                r.m  * ximesh.volLog(j.m ) + 

                r.ij * ximesh.volLog(j.ij) + 
                r.lj * ximesh.volLog(j.lj) + 
                r.lm * ximesh.volLog(j.lm) + 
                r.im * ximesh.volLog(j.im)
            );
}

template <typename XiMesh>
inline double interpVol(const Stencil<Cartesian, Wide, int>    j,
                        const Stencil<Cartesian, Wide, double> r,
                        const XiMesh& ximesh) {
    return  std::exp(
                r.o  * ximesh.volLog(j.o ) + 

                r.i  * ximesh.volLog(j.i ) + 
                r.j  * ximesh.volLog(j.j ) + 
                r.k  * ximesh.volLog(j.k ) + 
                r.l  * ximesh.volLog(j.l ) + 
                r.m  * ximesh.volLog(j.m ) + 
                r.n  * ximesh.volLog(j.n ) + 

                r.ij * ximesh.volLog(j.ij) + 
                r.lj * ximesh.volLog(j.lj) + 
                r.lm * ximesh.volLog(j.lm) + 
                r.im * ximesh.volLog(j.im) + 

                r.ik * ximesh.volLog(j.ik) + 
                r.lk * ximesh.volLog(j.lk) + 
                r.ln * ximesh.volLog(j.ln) + 
                r.in * ximesh.volLog(j.in) + 

                r.jk * ximesh.volLog(j.jk) + 
                r.mk * ximesh.volLog(j.mk) + 
                r.mn * ximesh.volLog(j.mn) + 
                r.jn * ximesh.volLog(j.jn)
            );
}

template <typename XiMesh>
inline double interpVol(const Stencil<Cylindrical, Tight, int>    j,
                        const Stencil<Cylindrical, Tight, double> r,
                        const XiMesh& ximesh) {
    return  std::exp(
                r.o  * ximesh.volLog(j.o ) + 
                r.x  * ximesh.volLog(j.x ) + 

                r.i  * ximesh.volLog(j.i ) + 
                r.j  * ximesh.volLog(j.j ) 
            );
}

template <typename XiMesh>
inline double interpVol(const Stencil<Cartesian, Tight, int>    j,
                        const Stencil<Cartesian, Tight, double> r,
                        const XiMesh& ximesh) {
    return  std::exp(
                r.o  * ximesh.volLog(j.o ) + 
                r.x  * ximesh.volLog(j.x ) + 

                r.i  * ximesh.volLog(j.i ) + 
                r.j  * ximesh.volLog(j.j ) + 
                r.k  * ximesh.volLog(j.k ) 
            );
}

template <typename XiMesh>
inline double interpVol(const Stencil<Cylindrical, Tight2, int>    j,
                        const Stencil<Cylindrical, Tight2, double> r,
                        const XiMesh& ximesh) {
    return  std::exp(
                r.o  * ximesh.volLog(j.o ) + 

                r.i  * ximesh.volLog(j.i ) + 
                r.j  * ximesh.volLog(j.j ) + 
                r.l  * ximesh.volLog(j.l ) + 
                r.m  * ximesh.volLog(j.m ) +

                r.ij * ximesh.volLog(j.ij)
            );
}

template <typename XiMesh>
inline double interpVol(const Stencil<Cartesian, Tight2, int>    j,
                        const Stencil<Cartesian, Tight2, double> r,
                        const XiMesh& ximesh) {
    return  std::exp(
                r.o  * ximesh.volLog(j.o ) + 

                r.i  * ximesh.volLog(j.i ) + 
                r.j  * ximesh.volLog(j.j ) + 
                r.k  * ximesh.volLog(j.k ) + 
                r.l  * ximesh.volLog(j.l ) + 
                r.m  * ximesh.volLog(j.m ) + 
                r.n  * ximesh.volLog(j.n ) +

                r.ij * ximesh.volLog(j.ij) +
                r.ik * ximesh.volLog(j.ik) +
                r.jk * ximesh.volLog(j.jk)
            );
}

template <typename XiMesh>
inline double interpVol(const Stencil<Cylindrical, Grad, int>    j,
                        const Stencil<Cylindrical, Grad, double> r,
                        const XiMesh& ximesh) {
    return  std::exp(
                r.o  * ximesh.volLog(j.o ) + 

                r.i  * ximesh.volLog(j.i ) + 
                r.j  * ximesh.volLog(j.j ) + 
                r.l  * ximesh.volLog(j.l ) + 
                r.m  * ximesh.volLog(j.m ) +

                r.ij * ximesh.volLog(j.ij) +

                r.ii * ximesh.volLog(j.ii) +
                r.jj * ximesh.volLog(j.jj)
            );
}

template <typename XiMesh>
inline double interpVol(const Stencil<Cartesian, Grad, int>    j,
                        const Stencil<Cartesian, Grad, double> r,
                        const XiMesh& ximesh) {
    return  std::exp(
                r.o  * ximesh.volLog(j.o ) + 

                r.i  * ximesh.volLog(j.i ) + 
                r.j  * ximesh.volLog(j.j ) + 
                r.k  * ximesh.volLog(j.k ) + 
                r.l  * ximesh.volLog(j.l ) + 
                r.m  * ximesh.volLog(j.m ) + 
                r.n  * ximesh.volLog(j.n ) +

                r.ij * ximesh.volLog(j.ij) +
                r.ik * ximesh.volLog(j.ik) +
                r.jk * ximesh.volLog(j.jk) +

                r.ii * ximesh.volLog(j.ii) +
                r.jj * ximesh.volLog(j.jj) +
                r.kk * ximesh.volLog(j.kk)
            );
}

inline const Stencil<Cylindrical, Symmetric, double> 
   operator*(Stencil<Cylindrical, Symmetric, double> y, const double x) {
    y.o  *= x;

    y.i  *= x;
    y.l  *= x;
    y.j  *= x;
    y.m  *= x;

    return y;
}

inline const Stencil<Cartesian, Symmetric, double> 
   operator*(Stencil<Cartesian, Symmetric, double> y, const double x) {
    y.o  *= x;

    y.i  *= x;
    y.l  *= x;
    y.j  *= x;
    y.m  *= x;
    y.k  *= x;
    y.n  *= x;

    return y;
}

inline const Stencil<Cylindrical, Wide, double> 
   operator*(Stencil<Cylindrical, Wide, double> y, const double x) {
    y.o  *= x;

    y.i  *= x;
    y.l  *= x;
    y.j  *= x;
    y.m  *= x;

    y.ij *= x;
    y.lj *= x;
    y.lm *= x;
    y.im *= x;

    return y;
}

inline const Stencil<Cartesian, Wide, double> 
   operator*(Stencil<Cartesian, Wide, double> y, const double x) {
    y.o  *= x;

    y.i  *= x;
    y.l  *= x;
    y.j  *= x;
    y.m  *= x;
    y.k  *= x;
    y.n  *= x;

    y.ij *= x;
    y.lj *= x;
    y.lm *= x;
    y.im *= x;

    y.ik *= x;
    y.lk *= x;
    y.ln *= x;
    y.in *= x;

    y.jk *= x;
    y.mk *= x;
    y.mn *= x;
    y.jn *= x;

    return y;
}

inline const Stencil<Cylindrical, Tight, double> 
   operator*(Stencil<Cylindrical, Tight, double> y, const double x) {
    y.o  *= x;
    y.x  *= x;

    y.i  *= x;
    y.j  *= x;

    return y;
}

inline const Stencil<Cartesian, Tight, double> 
   operator*(Stencil<Cartesian, Tight, double> y, const double x) {
    y.o  *= x;
    y.x  *= x;

    y.i  *= x;
    y.j  *= x;
    y.k  *= x;

    return y;
}

inline const Stencil<Cylindrical, Tight2, double> 
   operator*(Stencil<Cylindrical, Tight2, double> y, const double x) {
    y.o  *= x;

    y.i  *= x;
    y.l  *= x;
    y.j  *= x;
    y.m  *= x;

    y.ij *= x;

    return y;
}

inline const Stencil<Cartesian, Tight2, double> 
   operator*(Stencil<Cartesian, Tight2, double> y, const double x) {
    y.o  *= x;

    y.i  *= x;
    y.l  *= x;
    y.j  *= x;
    y.m  *= x;
    y.k  *= x;
    y.n  *= x;

    y.ij *= x;
    y.ik *= x;
    y.jk *= x;

    return y;
}

inline const Stencil<Cylindrical, Grad, double> 
   operator*(Stencil<Cylindrical, Grad, double> y, const double x) {
    y.o  *= x;

    y.i  *= x;
    y.l  *= x;
    y.j  *= x;
    y.m  *= x;

    y.ij *= x;

    y.ii *= x;
    y.jj *= x;

    return y;
}

inline const Stencil<Cartesian, Grad, double> 
   operator*(Stencil<Cartesian, Grad, double> y, const double x) {
    y.o  *= x;

    y.i  *= x;
    y.l  *= x;
    y.j  *= x;
    y.m  *= x;
    y.k  *= x;
    y.n  *= x;

    y.ij *= x;
    y.ik *= x;
    y.jk *= x;

    y.ii *= x;
    y.jj *= x;
    y.kk *= x;

    return y;
}

template <typename F>
inline void addF(F& q, Stencil<Cylindrical, Symmetric, int>    j,
                 const Stencil<Cylindrical, Symmetric, double> d) {
    q[j.o ] += d.o ;

    q[j.i ] += d.i ;
    q[j.l ] += d.l ;
    q[j.j ] += d.j ;
    q[j.m ] += d.m ;
}

template <typename F>
inline void addF(F& q, Stencil<Cartesian, Symmetric, int>    j,
                 const Stencil<Cartesian, Symmetric, double> d) {
    q[j.o ] += d.o ;

    q[j.i ] += d.i ;
    q[j.l ] += d.l ;
    q[j.j ] += d.j ;
    q[j.m ] += d.m ;
    q[j.k ] += d.k ;
    q[j.n ] += d.n ;
}

template <typename F>
inline void addF(F& q, Stencil<Cylindrical, Wide, int>    j,
                 const Stencil<Cylindrical, Wide, double> d) {
    q[j.o ] += d.o ;

    q[j.i ] += d.i ;
    q[j.l ] += d.l ;
    q[j.j ] += d.j ;
    q[j.m ] += d.m ;

    q[j.ij] += d.ij;
    q[j.lj] += d.lj;
    q[j.lm] += d.lm;
    q[j.im] += d.im;
}

template <typename F>
inline void addF(F& q, Stencil<Cartesian, Wide, int>    j,
                 const Stencil<Cartesian, Wide, double> d) {
    q[j.o ] += d.o ;

    q[j.i ] += d.i ;
    q[j.l ] += d.l ;
    q[j.j ] += d.j ;
    q[j.m ] += d.m ;
    q[j.k ] += d.k ;
    q[j.n ] += d.n ;

    q[j.ij] += d.ij;
    q[j.lj] += d.lj;
    q[j.lm] += d.lm;
    q[j.im] += d.im;
                   
    q[j.ik] += d.ik;
    q[j.lk] += d.lk;
    q[j.ln] += d.ln;
    q[j.in] += d.in;
                   
    q[j.jk] += d.jk;
    q[j.mk] += d.mk;
    q[j.mn] += d.mn;
    q[j.jn] += d.jn;
}

template <typename F>
inline void addF(F& q, Stencil<Cylindrical, Tight, int>    j,
                 const Stencil<Cylindrical, Tight, double> d) {
    q[j.o ] += d.o ;
    q[j.x ] += d.x ;

    q[j.i ] += d.i ;
    q[j.j ] += d.j ;
}

template <typename F>
inline void addF(F& q, Stencil<Cartesian, Tight, int>    j,
                 const Stencil<Cartesian, Tight, double> d) {
    q[j.o ] += d.o ;
    q[j.x ] += d.x ;

    q[j.i ] += d.i ;
    q[j.j ] += d.j ;
    q[j.k ] += d.k ;
}

template <typename F>
inline void addF(F& q, Stencil<Cylindrical, Tight2, int>    j,
                 const Stencil<Cylindrical, Tight2, double> d) {
    q[j.o ] += d.o ;

    q[j.i ] += d.i ;
    q[j.l ] += d.l ;
    q[j.j ] += d.j ;
    q[j.m ] += d.m ;

    q[j.ij] += d.ij;
}

template <typename F>
inline void addF(F& q, Stencil<Cartesian, Tight2, int>    j,
                 const Stencil<Cartesian, Tight2, double> d) {
    q[j.o ] += d.o ;

    q[j.i ] += d.i ;
    q[j.l ] += d.l ;
    q[j.j ] += d.j ;
    q[j.m ] += d.m ;
    q[j.k ] += d.k ;
    q[j.n ] += d.n ;

    q[j.ij] += d.ij;
    q[j.ik] += d.ik;
    q[j.jk] += d.jk;
}

template <typename F>
inline void addF(F& q, Stencil<Cylindrical, Grad, int>    j,
                 const Stencil<Cylindrical, Grad, double> d) {
    q[j.o ] += d.o ;

    q[j.i ] += d.i ;
    q[j.l ] += d.l ;
    q[j.j ] += d.j ;
    q[j.m ] += d.m ;

    q[j.ij] += d.ij;

    q[j.ii] += d.ii;
    q[j.jj] += d.jj;
}

template <typename F>
inline void addF(F& q, Stencil<Cartesian, Grad, int>    j,
                 const Stencil<Cartesian, Grad, double> d) {
    q[j.o ] += d.o ;

    q[j.i ] += d.i ;
    q[j.l ] += d.l ;
    q[j.j ] += d.j ;
    q[j.m ] += d.m ;
    q[j.k ] += d.k ;
    q[j.n ] += d.n ;

    q[j.ij] += d.ij;
    q[j.ik] += d.ik;
    q[j.jk] += d.jk;

    q[j.ii] += d.ii;
    q[j.jj] += d.jj;
    q[j.kk] += d.kk;
}

inline const Stencil<Cartesian, Symmetric, double>
         div(Stencil<Cartesian, Symmetric, double> x,
       const Stencil<Cartesian, Symmetric, double> y)
{
    x.o  /= y.o ;

    x.i  /= y.i ;
    x.l  /= y.l ;
    x.j  /= y.j ;
    x.m  /= y.m ;
    x.k  /= y.k ;
    x.n  /= y.n ;

    return x;
}

inline const Stencil<Cylindrical, Symmetric, double>
         div(Stencil<Cylindrical, Symmetric, double> x,
       const Stencil<Cylindrical, Symmetric, double> y)
{
    x.o  /= y.o ;

    x.i  /= y.i ;
    x.l  /= y.l ;
    x.j  /= y.j ;
    x.m  /= y.m ;

    return x;
}

inline const Stencil<Cartesian, Wide, double>
         div(Stencil<Cartesian, Wide, double> x,
       const Stencil<Cartesian, Wide, double> y)
{
    x.o  /= y.o ;

    x.i  /= y.i ;
    x.l  /= y.l ;
    x.j  /= y.j ;
    x.m  /= y.m ;
    x.k  /= y.k ;
    x.n  /= y.n ;

    x.ij /= y.ij;
    x.lj /= y.lj;
    x.lm /= y.lm;
    x.im /= y.im;

    x.ik /= y.ik;
    x.lk /= y.lk;
    x.ln /= y.ln;
    x.in /= y.in;

    x.jk /= y.jk;
    x.mk /= y.mk;
    x.mn /= y.mn;
    x.jn /= y.jn;

    return x;
}

inline const Stencil<Cylindrical, Wide, double>
         div(Stencil<Cylindrical, Wide, double> x,
       const Stencil<Cylindrical, Wide, double> y)
{
    x.o  /= y.o ;

    x.i  /= y.i ;
    x.l  /= y.l ;
    x.j  /= y.j ;
    x.m  /= y.m ;

    x.ij /= y.ij;
    x.lj /= y.lj;
    x.lm /= y.lm;
    x.im /= y.im;

    return x;
}

inline const Stencil<Cartesian, Tight, double>
         div(Stencil<Cartesian, Tight, double> x,
       const Stencil<Cartesian, Tight, double> y)
{
    x.o  /= y.o ;
    x.x  /= y.x ;

    x.i  /= y.i ;
    x.j  /= y.j ;
    x.k  /= y.k ;

    return x;
}

inline const Stencil<Cylindrical, Tight, double>
         div(Stencil<Cylindrical, Tight, double> x,
       const Stencil<Cylindrical, Tight, double> y)
{
    x.o  /= y.o ;
    x.x  /= y.x ;

    x.i  /= y.i ;
    x.j  /= y.j ;

    return x;
}

inline const Stencil<Cartesian, Tight2, double>
         div(Stencil<Cartesian, Tight2, double> x,
       const Stencil<Cartesian, Tight2, double> y)
{
    x.o  /= y.o ;

    x.i  /= y.i ;
    x.l  /= y.l ;
    x.j  /= y.j ;
    x.m  /= y.m ;
    x.k  /= y.k ;
    x.n  /= y.n ;

    x.ij /= y.ij;
    x.ik /= y.ik;
    x.jk /= y.jk;

    return x;
}

inline const Stencil<Cylindrical, Tight2, double>
         div(Stencil<Cylindrical, Tight2, double> x,
       const Stencil<Cylindrical, Tight2, double> y)
{
    x.o  /= y.o ;

    x.i  /= y.i ;
    x.l  /= y.l ;
    x.j  /= y.j ;
    x.m  /= y.m ;

    x.ij /= y.ij;

    return x;
}

inline const Stencil<Cartesian, Grad, double>
         div(Stencil<Cartesian, Grad, double> x,
       const Stencil<Cartesian, Grad, double> y)
{
    x.o  /= y.o ;

    x.i  /= y.i ;
    x.l  /= y.l ;
    x.j  /= y.j ;
    x.m  /= y.m ;
    x.k  /= y.k ;
    x.n  /= y.n ;

    x.ij /= y.ij;
    x.ik /= y.ik;
    x.jk /= y.jk;

    x.ii /= y.ii;
    x.jj /= y.jj;
    x.kk /= y.kk;

    return x;
}

inline const Stencil<Cylindrical, Grad, double>
         div(Stencil<Cylindrical, Grad, double> x,
       const Stencil<Cylindrical, Grad, double> y)
{
    x.o  /= y.o ;

    x.i  /= y.i ;
    x.l  /= y.l ;
    x.j  /= y.j ;
    x.m  /= y.m ;

    x.ij /= y.ij;

    x.ii /= y.ii;
    x.jj /= y.jj;

    return x;
}

template <PowMethod powmethod> inline
void log_(const double* p1, double* p2);

template <> inline
void log_<Std>(const double* p1, double* p2) {
    *p2       = std::log(*p1);
    *(p2 + 1) = std::log(*(p1 + 1));
}

template <> inline
void log_<FastSSE>(const double* p1, double* p2) {
    __m128d md;
    md = _mm_load_pd(p1);
    md = sse::log(md);
    _mm_store_pd(p2, md);
}

template <PowMethod powmethod>
inline const Stencil<Cartesian, Symmetric, double>
       mylog(Stencil<Cartesian, Symmetric, double> x)
{
    log_<powmethod>(&x.o, &x.o);
    log_<powmethod>(&x.i, &x.i);
    log_<powmethod>(&x.j, &x.j);
    log_<powmethod>(&x.k, &x.k);
    return x;
}

template <PowMethod powmethod>
inline const Stencil<Cylindrical, Symmetric, double>
       mylog(Stencil<Cylindrical, Symmetric, double> x)
{
    log_<powmethod>(&x.o, &x.o);
    log_<powmethod>(&x.i, &x.i);
    log_<powmethod>(&x.j, &x.j);
    return x;
}

template <PowMethod powmethod>
inline const Stencil<Cartesian, Wide, double>
       mylog(Stencil<Cartesian, Wide, double> x)
{
    log_<powmethod>(&x.o, &x.o);
    log_<powmethod>(&x.i, &x.i);
    log_<powmethod>(&x.j, &x.j);
    log_<powmethod>(&x.k, &x.k);

    log_<powmethod>(&x.ij, &x.ij);
    log_<powmethod>(&x.im, &x.im);
                                
    log_<powmethod>(&x.ik, &x.ik);
    log_<powmethod>(&x.in, &x.in);
                                
    log_<powmethod>(&x.jk, &x.jk);
    log_<powmethod>(&x.jn, &x.jn);
    return x;
}

template <PowMethod powmethod>
inline const Stencil<Cylindrical, Wide, double>
       mylog(Stencil<Cylindrical, Wide, double> x)
{
    log_<powmethod>(&x.o, &x.o);
    log_<powmethod>(&x.i, &x.i);
    log_<powmethod>(&x.j, &x.j);

    log_<powmethod>(&x.ij, &x.ij);
    log_<powmethod>(&x.im, &x.im);
    return x;
}

template <PowMethod powmethod>
inline const Stencil<Cartesian, Tight, double>
       mylog(Stencil<Cartesian, Tight, double> x)
{
    log_<powmethod>(&x.o, &x.o);
    log_<powmethod>(&x.i, &x.i);
    log_<powmethod>(&x.k, &x.k);
    return x;
}

template <PowMethod powmethod>
inline const Stencil<Cylindrical, Tight, double>
       mylog(Stencil<Cylindrical, Tight, double> x)
{
    log_<powmethod>(&x.o, &x.o);
    log_<powmethod>(&x.i, &x.i);
    return x;
}

template <PowMethod powmethod>
inline const Stencil<Cartesian, Tight2, double>
       mylog(Stencil<Cartesian, Tight2, double> x)
{
    log_<powmethod>(&x.i,  &x.i);
    log_<powmethod>(&x.j,  &x.j);
    log_<powmethod>(&x.k,  &x.k);
    log_<powmethod>(&x.o,  &x.o);
    log_<powmethod>(&x.ik, &x.ik);
    return x;
}

template <PowMethod powmethod>
inline const Stencil<Cylindrical, Tight2, double>
       mylog(Stencil<Cylindrical, Tight2, double> x)
{
    log_<powmethod>(&x.o, &x.o);
    log_<powmethod>(&x.i, &x.i);
    log_<powmethod>(&x.j, &x.j);
    return x;
}

template <PowMethod powmethod>
inline const Stencil<Cartesian, Grad, double>
       mylog(Stencil<Cartesian, Grad, double> x)
{
    log_<powmethod>(&x.o,  &x.o);
    log_<powmethod>(&x.i,  &x.i);
    log_<powmethod>(&x.j,  &x.j);
    log_<powmethod>(&x.k,  &x.k);
    log_<powmethod>(&x.ij, &x.ij);
    log_<powmethod>(&x.jk, &x.jk);
    log_<powmethod>(&x.jj, &x.jj);
    return x;
}

template <PowMethod powmethod>
inline const Stencil<Cylindrical, Grad, double>
       mylog(Stencil<Cylindrical, Grad, double> x)
{
    log_<powmethod>(&x.o,  &x.o);
    log_<powmethod>(&x.i,  &x.i);
    log_<powmethod>(&x.j,  &x.j);
    log_<powmethod>(&x.ii, &x.ii);
    return x;
}

template <Volume volume, Symmetry symmetry, template <Symmetry> class XiMeshType>
struct MakeRSwarm;

template <template <Symmetry> class XiMeshType>
struct MakeRSwarm<Symmetric, Cylindrical, XiMeshType> {
    static const Symmetry symmetry = Cylindrical;
    static const Volume volume = Symmetric;
    const Stencil<symmetry, volume, double> 
            make(const typename SymmetryTrait<symmetry>::Vd x,
                 const Stencil<symmetry, volume, int> j,
                 const XiMeshType<symmetry>& ximesh) const
    {
        Stencil<symmetry, volume, double> s;

//        std::cout << x << ' ' << t.b << ' ' << t.f << std::endl;

        V2d hf, hb;

        hf[0] =   0.5 * ( ximesh.vvoli(j.i, 0) + ximesh.vvoli(j.o, 0) );
        hb[0] = - 0.5 * ( ximesh.vvoli(j.l, 0) + ximesh.vvoli(j.o, 0) );
        hf[1] =   0.5 * ( ximesh.vvoli(j.j, 1) + ximesh.vvoli(j.o, 1) );
        hb[1] = - 0.5 * ( ximesh.vvoli(j.m, 1) + ximesh.vvoli(j.o, 1) );
        
/*
        std::cout << ximesh[j.i][0] - ximesh[j.o][0] << " = " << hf[0] << std::endl;
        std::cout << ximesh[j.l][0] - ximesh[j.o][0] << " = " << hb[0] << std::endl;
        std::cout << ximesh[j.j][1] - ximesh[j.o][1] << " = " << hf[1] << std::endl;
        std::cout << ximesh[j.m][1] - ximesh[j.o][1] << " = " << hb[1] << std::endl;
*/

        s.i  =   x[0] * (x[0] - hb[0]) / hf[0] / (hf[0] - hb[0]);
        s.l  = - x[0] * (x[0] - hf[0]) / hb[0] / (hf[0] - hb[0]);
        s.j  =   x[1] * (x[1] - hb[1]) / hf[1] / (hf[1] - hb[1]);
        s.m  = - x[1] * (x[1] - hf[1]) / hb[1] / (hf[1] - hb[1]);

        s.o  =  1 - s.i - s.j - s.l - s.m;

//        std::cout << s << std::endl;

        return s;
    }
};

template <template <Symmetry> class XiMeshType>
struct MakeRSwarm<Symmetric, Cartesian, XiMeshType> {
    static const Symmetry symmetry = Cartesian;
    static const Volume volume = Symmetric;
    const Stencil<symmetry, volume, double> 
            make(const typename SymmetryTrait<symmetry>::Vd x,
                 const Stencil<symmetry, volume, int> j,
                 const XiMeshType<symmetry>& ximesh) const
    {
        Stencil<symmetry, volume, double> s;

        V3d hf, hb;

        hf[0] =   0.5 * ( ximesh.vvoli(j.i, 0) + ximesh.vvoli(j.o, 0) );
        hb[0] = - 0.5 * ( ximesh.vvoli(j.l, 0) + ximesh.vvoli(j.o, 0) );
        hf[1] =   0.5 * ( ximesh.vvoli(j.j, 1) + ximesh.vvoli(j.o, 1) );
        hb[1] = - 0.5 * ( ximesh.vvoli(j.m, 1) + ximesh.vvoli(j.o, 1) );
        hf[2] =   0.5 * ( ximesh.vvoli(j.k, 2) + ximesh.vvoli(j.o, 2) );
        hb[2] = - 0.5 * ( ximesh.vvoli(j.n, 2) + ximesh.vvoli(j.o, 2) );

        s.i  =   x[0] * (x[0] - hb[0]) / hf[0] / (hf[0] - hb[0]);
        s.l  = - x[0] * (x[0] - hf[0]) / hb[0] / (hf[0] - hb[0]);
        s.j  =   x[1] * (x[1] - hb[1]) / hf[1] / (hf[1] - hb[1]);
        s.m  = - x[1] * (x[1] - hf[1]) / hb[1] / (hf[1] - hb[1]);
        s.k  =   x[2] * (x[2] - hb[2]) / hf[2] / (hf[2] - hb[2]);
        s.n  = - x[2] * (x[2] - hf[2]) / hb[2] / (hf[2] - hb[2]);

        s.o  =  1 - s.i - s.j - s.k - s.l - s.m - s.n;

        return s;
    }
};
template <template <Symmetry> class XiMeshType>
struct MakeRSwarm<Tight, Cylindrical, XiMeshType> {
    static const Symmetry symmetry = Cylindrical;
    static const Volume volume = Tight;
    const Stencil<symmetry, volume, double> 
            make(const typename SymmetryTrait<symmetry>::Vd x,
                 const Stencil<symmetry, volume, int> j,
                 const XiMeshType<symmetry>& ximesh) const
    {
        Stencil<symmetry, volume, double> s;

        V2d hf, hb;

        hf[0] =   0.5 * ( ximesh.vvoli(j.i, 0) + ximesh.vvoli(j.o, 0) );
        hb[0] = - 0.5 * ( ximesh.vvoli(j.x, 0) + ximesh.vvoli(j.o, 0) );
        hf[1] =   0.5 * ( ximesh.vvoli(j.j, 1) + ximesh.vvoli(j.o, 1) );
        hb[1] = - 0.5 * ( ximesh.vvoli(j.x, 1) + ximesh.vvoli(j.o, 1) );

        V2d ax = abs(x);

        s.x = dot(ax, hf - ax) / dot(hb, hf - hb);

        V2d sij = (ax - hb * s.x) / hf;
        
        s.i  =  sij[0];
        s.j  =  sij[1];

        s.o  =  1 - s.i - s.j - s.x;

        return s;
    }
};

template <template <Symmetry> class XiMeshType>
struct MakeRSwarm<Tight, Cartesian, XiMeshType> {
    static const Symmetry symmetry = Cartesian;
    static const Volume volume = Tight;
    const Stencil<symmetry, volume, double> 
            make(const typename SymmetryTrait<symmetry>::Vd x,
                 const Stencil<symmetry, volume, int> j,
                 const XiMeshType<symmetry>& ximesh) const
    {
        Stencil<symmetry, volume, double> s;

        V3d hf, hb;

        hf[0] =   0.5 * ( ximesh.vvoli(j.i, 0) + ximesh.vvoli(j.o, 0) );
        hb[0] = - 0.5 * ( ximesh.vvoli(j.x, 0) + ximesh.vvoli(j.o, 0) );
        hf[1] =   0.5 * ( ximesh.vvoli(j.j, 1) + ximesh.vvoli(j.o, 1) );
        hb[1] = - 0.5 * ( ximesh.vvoli(j.x, 1) + ximesh.vvoli(j.o, 1) );
        hf[2] =   0.5 * ( ximesh.vvoli(j.k, 2) + ximesh.vvoli(j.o, 2) );
        hb[2] = - 0.5 * ( ximesh.vvoli(j.x, 2) + ximesh.vvoli(j.o, 2) );

        V3d ax = abs(x);

        s.x = dot(ax, hf - ax) / dot(hb, hf - hb);

        V3d sijk = (ax - hb * s.x) / hf;
        
        s.i  =  sijk[0];
        s.j  =  sijk[1];
        s.k  =  sijk[2];

        s.o  =  1 - s.i - s.j - s.k - s.x;

        return s;
    }
};


template <Volume volume, Symmetry symmetry, template <Symmetry> class XiMeshType>
Stencil<symmetry, volume, double> doMakeRSwarm(const typename SymmetryTrait<symmetry>::Vd x,
                                             const Stencil<symmetry, volume, int> j,
                                             const XiMeshType<symmetry>& ximesh)
{
    return MakeRSwarm<volume, symmetry, XiMeshType>().make(x, j, ximesh);
}


#endif
