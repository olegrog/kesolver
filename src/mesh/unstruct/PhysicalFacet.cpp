#include "PhysicalFacet.hpp"
#include "math.h"

inline V3d triangleCenter(V3d p1, V3d p2, V3d p3) {
    return (p1 + p2 + p3) / 3.;
}

inline double triangleVolume(V3d p1, V3d p2, V3d p3) {
    return 0.5 * norm( cross(p2 - p1, p3 - p1) );
}

void PhysicalFacet::findNormalAndSquare()
{
    n = cross(vertex[1] - vertex[0], vertex[2] - vertex[0]);
    n /= norm(n);
    if (type == 2) {
        S = triangleVolume(vertex[0], vertex[1], vertex[2]);
        center = triangleCenter(vertex[0], vertex[1], vertex[2]);
    }
    else if (type == 3) {
        double s1 =     triangleVolume(vertex[0], vertex[1], vertex[2]);
        V3d c1 =        triangleCenter(vertex[0], vertex[1], vertex[2]);
        double s2 =     triangleVolume(vertex[0], vertex[2], vertex[3]);
        V3d c2 =        triangleCenter(vertex[0], vertex[2], vertex[3]);
        S = s1 + s2;
        center = (c1 * s1 + c2 * s2) / (s1 + s2);
    }
}

void PhysicalFacet::setType(int type_)
{
    type = type_;
    if(type == 2) numberOfVertex = 3;
    if(type == 3) numberOfVertex = 4;
}

void PhysicalFacet::setNeigbors(const std::vector<int>& neigbors)
{
    polygon = neigbors;
}

void PhysicalFacet::setVertexes(const std::vector<V3d>& vertexes)
{
    vertex = vertexes;
    findNormalAndSquare();
}

void PhysicalFacet::findMultInOut(double t, const std::vector<Polygon*>& spacemesh) {
    dt = t;
    mult_in = t * S / spacemesh[polygon[0]]->getVolume();
    if (polygon.size() > 1)
        mult_out = t * S / spacemesh[polygon[1]]->getVolume();
}

void PhysicalFacet::findPhi(const std::vector<Polygon*>& spacemesh, const Gas& gas)
{
    doFindPhi(spacemesh, gas);
}
void PhysicalFacet::transfer(std::vector<Polygon*>& spacemesh,
        const Gas& gas) 
{
/*
    for (size_t j = 0; j < polygon.size(); ++j) {
        double f = spacemesh[polygon[j]]->f().g()[0];
        if (f != f) {
            LABEL
            std::cout << polygon.size() << std::endl;
            std::cout << "nan bef " << polygon[j] << std::endl;
            exit(-1);
        }
    }
*/
    doTransfer(spacemesh, gas);
/*
    for (size_t j = 0; j < polygon.size(); ++j) {
        double f = spacemesh[polygon[j]]->f().g()[0];
        if (f != f) {
            LABEL
            std::cout << polygon.size() << std::endl;
            std::cout << "nan aft " << polygon[j] << std::endl;
            exit(-1);
        }
    }
*/
}

void PhysicalFacet::transfer2(std::vector<Polygon*>& spacemesh,
                              const Gas& gas) 
{
    doTransfer2(spacemesh, gas);
}
void PhysicalFacet::findGradient(const std::vector<Polygon*>& spacemesh,
                                 const Gas& gas) 
{
    doFindGradient(spacemesh, gas);
}

void PhysicalFacet::doFindPhi(const std::vector<Polygon*>& spacemesh,
                              const Gas& gas)
{
    for (size_t j = 0; j < polygon.size(); ++j) {
        Polygon* poly = spacemesh[polygon[j]];
        SpeedFunction& function = poly->f();
        V3d cc = getCenter() - poly->getCenter();
        for(size_t i = 0; i < function.size(); i++) {
            double f = function.f()[i];
            V3d df = function.getGradient()[i];
            double dd = - 0.5 * dt * gas.dot(i, df);
//            double dd = 0.0;
            double dl = dot(df, cc + dd);
            double d1 = function.getFMax()[i] - f;
            double d2 = function.getFMin()[i] - f;
            double phi;
            if (std::abs(dl) < 1e-10) phi = 0.;
            else {
                if (dl > 0) phi = std::min(1., d1 / dl);
                else        phi = std::min(1., d2 / dl);
            }

            if (function.getPhi()[i] > phi)
                function.getPhi()[i] = phi;
            
        }
    }
}

