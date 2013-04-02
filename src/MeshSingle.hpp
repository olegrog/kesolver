#ifndef _MESHSINGLE_HPP_
#define _MESHSINGLE_HPP_

#include <vector>

#include "Mesh.hpp"

class MeshSingle : public Mesh {
    public:
        MeshSingle(int argc, char** argv, MeshBase* mesh);
            
        Cells&  getAllCells()      { return cells; }
        Cells&  getFlowingCells()  { return cells; }
        Cells&  getMyCells()       { return cells; }

        Facets& getAllFacets()     { return facets; }
        Facets& getFlowingFacets() { return facets; }

        const Ints& getMyCellIndexes() const { return my_cell_indexes; }

        double getTimeStep() const {
            return time_step;
        //    mesh_ptr->getTimeStep();
        }

        void newStep() {}

        int getRank() const { return 0; }
        int getSize() const { return 1; }

        void barrier() const {}

    private:
        double time_step;

        Mesh*   mesh_ptr;  

        Cells&  cells;
        Ints    my_cell_indexes;

        Facets& facets;
};

#endif // _MESHSINGLE_HPP_

