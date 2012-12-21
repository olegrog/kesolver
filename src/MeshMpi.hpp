#ifndef _MESHMPI_HPP_
#define _MESHMPI_HPP_

#include <vector>

#include "mesh/Mesh.hpp"

class MeshMpi {
    public:
        MeshMpi(const Mesh* mesh_ptr, int rank);

        typedef std::vector<int> Ints;

        Ints&  getMyCells()  { return mycells;  }
        Ints&  getMyFacets() { return myfacets; }

        // TODO: add all neccessary methods

    private:
        Ints   mycells;  
        Ints   myfacets;  
};

#endif // _MESHMPI_HPP_

