#include "mpi.h"

#include "MeshMpi.hpp"

MeshMpi::MeshMpi(MeshBase* mesh) :
    time_step(mesh->getTimeStep()),
    mpi_init(false),
    cells(mesh->getCells()),
    facets(mesh->getFacets())
{
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::vector<bool> is_cell_flowing(cells.size(), false);

    for (size_t i = 0; i < cells.size(); ++i) 
        if (cells[i]->getRank() == rank) {
            my_cells.push_back(cells[i]);
            my_cell_indexes.push_back(i);
            flowing_cells.push_back(cells[i]);
            is_cell_flowing[i] = true;
        }
    
    data_exchanger.init2(cells, rank); // TODO

	const std::vector<int>& exchange_cells =
            data_exchanger.getExchangeCells();

    for (size_t i = 0; i < exchange_cells.size(); ++i) {
        int j = exchange_cells[i];
        flowing_cells.push_back(cells[j]);
        is_cell_flowing[j] = true;
    }

    for (size_t i = 0; i < facets.size(); ++i) {
        bool push = true;
        for (size_t j = 0; j < facets[i]->getNeigbors().size(); ++j) 
            push &= is_cell_flowing[facets[i]->getNeigbors()[j]];
        if (push)
            flowing_facets.push_back(facets[i]);
    }

}

void MeshMpi::newStep()
{
    if (not mpi_init)
        data_exchanger.mpiInit(cells); 
	data_exchanger.swap();
}

void MeshMpi::barrier() const 
{
    MPI_Barrier(MPI_COMM_WORLD);
}

MeshMpi::~MeshMpi()
{
    MPI_Finalize();
}

