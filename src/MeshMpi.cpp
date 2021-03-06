#include "mpi.h"

#include "MeshMpi.hpp"

MeshMpi::MeshMpi(int argc, char** argv, MeshBase* mesh) :
    time_step(mesh->getTimeStep()),
    mpi_init(false),
    cells(mesh->getCells()),
    facets(mesh->getFacets())
{
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::vector<bool> is_cell_flowing(cells.size(), false);
    std::vector<bool> is_cell_my(cells.size(), false);

    for (size_t i = 0; i < cells.size(); ++i) 
        if (cells[i]->getRank() == rank) {
            my_cells.push_back(cells[i]);
            my_cell_indexes.push_back(i);
            flowing_cells.push_back(cells[i]);
            is_cell_flowing[i] = true;
            is_cell_my[i] = true;
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
        bool is_flowing = true;
        for (size_t j = 0; j < facets[i]->getNeigbors().size(); ++j) 
            is_flowing &= is_cell_flowing[facets[i]->getNeigbors()[j]];
        if (is_flowing)
            flowing_facets.push_back(facets[i]);
        bool is_my = false;
        for (size_t j = 0; j < facets[i]->getNeigbors().size(); ++j) 
            is_my |= is_cell_my[facets[i]->getNeigbors()[j]];
        if (is_my)
            my_facets.push_back(facets[i]);
    }

}

void MeshMpi::newStep()
{
    if (not mpi_init) {
        data_exchanger.mpiInit(cells); 
        mpi_init = true;
    }
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

