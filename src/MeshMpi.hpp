#ifndef _MESHMPI_HPP_
#define _MESHMPI_HPP_

#include <vector>

#include "mesh/Mesh.hpp"
#include "DataExchange.hpp"

class MeshMpi {
    public:
        MeshMpi(const Mesh* mesh_ptr);

        typedef typename Mesh::Cells  Cells;
        typedef typename Mesh::Facets Facets;

        Cells&  getAllCells()      { return cells; }
        Cells&  getFlowingCells()  { return flowing_cells; }
        Cells&  getMyCells()       { return my_cells; }

        Facets& getAllFacets()     { return facets; }
        Facets& getFlowingFacets() { return flowing_facets; }

        void newStep();

    private:
        bool mpi_init;

        Mesh*   mesh_ptr;  

        Cells&  cells;
        Cells   flowing_cells;  
        Cells   my_cells;  

        Facets& facets;
        Facets  my_facets;  

		DataExchanger data_exchanger;

};

#endif // _MESHMPI_HPP_

