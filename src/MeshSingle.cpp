#include "MeshSingle.hpp"

MeshSingle::MeshSingle(MeshBase* mesh) :
    time_step(mesh->getTimeStep()),
    cells(mesh->getCells()),
    facets(mesh->getFacets())
{
    my_cell_indexes.reserve(cells.size());
    for (size_t i = 0; i < cells.size(); ++i) 
        my_cell_indexes.push_back(i);
}
