#ifndef _UNSTRUCTMESH_HPP_
#define _UNSTRUCTMESH_HPP_

#include <vector>

#include "Polygon.hpp"
#include "PhysicalFacet.hpp"

#include "property_tree/property_tree.hpp"

class UnstructMesh {
    public:
        UnstructMesh(const PropertyTree& tree, const Gas& gas);

        void setMpiRank(int rank);

        typedef std::vector<Polygon*>       Cells;
        typedef std::vector<PhysicalFacet*> Facets;

        Cells&  getCells()  { return cells;  }
        Facets& getFacets() { return facets; }

        typedef std::vector<int> Ints;

        Ints&  getMyCells()  { return mycells;  }
        Ints&  getMyFacets() { return myfacets; }

        double getTimeStep() const { return time_step; }

    private:
        Cells  cells;                   
        Facets facets;  
        Ints   mycells;  
        Ints   myfacets;  

        double time_step;

};

#endif // _UNSTRUCTMESH_HPP_
