#ifndef _MESHMPI_HPP_
#define _MESHMPI_HPP_

#include <vector>

#include "Mesh.hpp"

#include "DataExchanger.hpp"

class MeshMpi : public Mesh {
    public:
        MeshMpi(int argc, char** argv, MeshBase* mesh_ptr);
        ~MeshMpi();

        Cells&  getAllCells()      { return cells; }
        Cells&  getFlowingCells()  { return flowing_cells; }
        Cells&  getMyCells()       { return my_cells; }

        Facets& getAllFacets()     { return facets; }
        Facets& getFlowingFacets() { return flowing_facets; }

        const Ints& getMyCellIndexes() const { return my_cell_indexes; }

        double getTimeStep() const {
            return time_step;
        //    mesh_ptr->getTimeStep();
        }

        void newStep();

        int getRank() const { return rank; }
        int getSize() const { return size; }

        void barrier() const;

    private:
        double time_step;
        bool mpi_init;

        Mesh*   mesh_ptr;  

        Cells&  cells;
        Cells   flowing_cells;  
        Cells   my_cells;  
        Ints    my_cell_indexes;

        Facets& facets;
        Facets  flowing_facets;  

		DataExchanger data_exchanger;
        int rank, size;

};

#endif // _MESHMPI_HPP_

