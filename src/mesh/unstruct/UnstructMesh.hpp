#ifndef _UNSTRUCTMESH_HPP_
#define _UNSTRUCTMESH_HPP_

#include <vector>

#include "mesh/Mesh.hpp"

#include "property_tree/property_tree.hpp"

class UnstructMesh : public Mesh {
    public:
        UnstructMesh(const PropertyTree& tree, const Gas& gas);

        virtual Cells&  getCells()  { return cells;  }
        virtual Facets& getFacets() { return facets; }

        virtual double getTimeStep() const {
            LABEL
            return time_step;
        }

    private:
        Cells  cells;                   
        Facets facets;  

        double time_step;
};

#endif // _UNSTRUCTMESH_HPP_
