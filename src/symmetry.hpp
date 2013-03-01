#ifndef _SYMMETRY_HPP_
#define _SYMMETRY_HPP_

#include "v.hpp"

enum Symmetry {
    Cylindrical,
    Cartesian
};

inline std::ostream& operator<<(std::ostream& to, Symmetry symmetry) {
    if       (symmetry == Cylindrical) return to << "Cylindrical";
    else if  (symmetry == Cartesian)   return to << "Cartesian"  ;
    else                               return to << symmetry     ;
}

template <Symmetry> struct SymmetryTrait;

template <>
struct SymmetryTrait<Cartesian> {
    typedef V3d Vd;
    typedef V3d Vm;
    typedef V3i Vi;
    typedef V3i Vj;
    typedef V3< std::vector<double> > VVd;
    static const Vd t_mult() { return Vd(1., 1., 1.); }
    static const int ximesh_dim = 3;
};

template <>
struct SymmetryTrait<Cylindrical> {
    typedef V2d    Vd;
    typedef double Vm;
    typedef V2i    Vi;
    typedef int    Vj;
    typedef V2< std::vector<double> > VVd;
    static const Vd t_mult() { return Vd(1., 1./2.); }
    static const int ximesh_dim = 2;
};

inline SymmetryTrait<Cartesian>::Vd vm2vd(SymmetryTrait<Cartesian>::Vm x) {
    return x;
}
inline SymmetryTrait<Cylindrical>::Vd vm2vd(SymmetryTrait<Cylindrical>::Vm x) {
    return SymmetryTrait<Cylindrical>::Vd(x, 0.0);
}

template <Symmetry symmetry>
inline typename SymmetryTrait<symmetry>::Vm x2vm(const double x);

template <>
inline SymmetryTrait<Cartesian>::Vm x2vm<Cartesian>(const double x) {
    return V3d(x, 0, 0);
}

template <>
inline SymmetryTrait<Cylindrical>::Vm x2vm<Cylindrical>(const double x) {
    return x;
}

#endif
