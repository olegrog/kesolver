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

void PhysicalFacet::setPolygonNumbers(const std::vector<int>& neigbors_)
{
    polygon = neigbors_;
}

void PhysicalFacet::setVertex(const std::vector<V3d>& vertex_)
{
    vertex = vertex_;
    findNormalAndSquare();
}

void PhysicalFacet::findMultInOut(double t, const std::vector<Polygon*>& spacemesh) {
    mult_in = t * S / spacemesh[polygon[0]]->getVolume();
    if (polygon.size() > 1)
        mult_out = t * S / spacemesh[polygon[1]]->getVolume();
}

void PhysicalFacet::findPhi(const std::vector<Polygon*>& spacemesh)
{
    if (is_active)
        doFindPhi(spacemesh);
}
void PhysicalFacet::transfer(std::vector<Polygon*>& spacemesh,
        const Gas& gas) 
{
    if (is_active) {
        for (size_t j = 0; j < polygon.size(); ++j) {
            double f = spacemesh[polygon[j]]->f().g()[0];
            if (f != f) {
                LABEL
                std::cout << polygon.size() << std::endl;
                std::cout << "nan bef " << polygon[j] << std::endl;
                exit(-1);
            }
        }

        doTransfer(spacemesh, gas);

        for (size_t j = 0; j < polygon.size(); ++j) {
            double f = spacemesh[polygon[j]]->f().g()[0];
            if (f != f) {
                LABEL
                std::cout << polygon.size() << std::endl;
                std::cout << "nan aft " << polygon[j] << std::endl;
                exit(-1);
            }
        }
    }
}
void PhysicalFacet::transfer2(std::vector<Polygon*>& spacemesh,
        const Gas& gas) 
{
    if (is_active)
        doTransfer2(spacemesh, gas);
}
void PhysicalFacet::findGradient(const std::vector<Polygon*>& spacemesh,
        const Gas& gas) 
{
    if (is_active)
        doFindGradient(spacemesh, gas);
}

void PhysicalFacet::doFindPhi(const std::vector<Polygon*>& spacemesh)
{
    for (size_t j = 0; j < polygon.size(); ++j) {
        Polygon* poly = spacemesh[polygon[j]];
        SpeedFunction& function = poly->f();
        for(size_t i = 0; i < function.size(); i++) {
            double f = function.f()[i];
            double p = f + dot(function.getGradient()[i], 
                    getCenter() - poly->getCenter());
            double d1 = function.getFMax()[i] - f;
            double d2 = function.getFMin()[i] - f;
            double dl = p - f;
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

void PhysicalFacet::activation(const std::vector<Polygon*>& spacemesh) {
    is_active = true;
    for (size_t i = 0; i < polygon.size(); ++i) 
        is_active = is_active && spacemesh[polygon[i]]->isActive();
}

