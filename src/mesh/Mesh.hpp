#ifndef _MESH_HPP_
#define _MESH_HPP_

#include <vector>

#include "unstruct/Polygon.hpp"
#include "unstruct/PhysicalFacet.hpp"

class Mesh {
    public:
        typedef std::vector<Polygon*>       Cells;
        typedef std::vector<PhysicalFacet*> Facets;

        virtual Cells&  getCells()  = 0;
        virtual Facets& getFacets() = 0;

        virtual double getTimeStep() const  = 0;
};

#endif // _MESH_HPP_

