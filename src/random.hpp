#ifndef _RANDOM_HPP_
#define _RANDOM_HPP_

#include "v.hpp"
#include "auxiliary.hpp"

inline const V3d randomSphere(const V2d point) {
    double phi = point[0] * 2 * M_PI;
    double ct  = 2 * point[1] - 1;
    double st  = std::sqrt( 1 - sqr(ct) );
    return V3d( st * std::cos(phi), st * std::sin(phi), ct );
}

inline const V2d randomHalfCircle(const V2d point) {
    double phi = point[0] * M_PI;
    double r  = std::sqrt(point[1]);
    return V2d( r * std::cos(phi), r * std::sin(phi) );
}

inline const V3d randomBall(const V3d point) {
    double r = cbrt( point[0] );
    return randomSphere(V2d(point[1], point[2])) * r;
}

#endif
