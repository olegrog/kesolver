#include <iostream>
#include <vector>
#include <string>

#include "mpi.h"

#include "property_tree/property_tree.hpp"
#include "logger/logger.hpp"

#include "Constructors.hpp"

#include "mesh/unstruct/UnstructMesh.hpp"

#include "Transfer.hpp"
#include "Transfer1.hpp"
#include "Transfer2.hpp"

#include "Integral.hpp"

#include "Printer.hpp"

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    PropertyTree prop_tree(argv[1]);

    Gas* gas_p;
    GasConstructor(prop_tree["gas"], &gas_p);
    Gas& gas = *gas_p;

    LOG(INFO) << "gas.size() = " << gas_p->size();

    UnstructMesh mesh(prop_tree, gas);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    mesh.setMpiRank(rank);

    LOG(INFO) << "rank = " << rank 
              << " facets.size() = "   << mesh.getFacets().size() 
              << " mypolys.size() = "  << mesh.getMyCells().size()
              << " myfacets.size() = " << mesh.getMyFacets().size();

    Transfer* transfer;
    int order = prop_tree["transfer"].isMember("order") ? 
                prop_tree["transfer"]["order"].asInt() : 1;
    LOG(INFO) << "order = " << order;

    if (order == 2)
        transfer = new Transfer2(mesh.getFacets(), mesh.getCells(), mesh.getMyCells());
    else
        transfer = new Transfer1();

    transfer->init(prop_tree, gas, 
                   mesh.getFacets(), mesh.getCells(), mesh.getMyCells(),
                   rank);

    Integral integral = IntegralConstructor(prop_tree["integral"]);

    int rep = prop_tree["transfer"].isMember("rep") ? 
              prop_tree["transfer"]["rep"].asInt() : 1;
    LOG(INFO) << "rep = " << rep;

    Printer printer(prop_tree["printer"]);
    for (int i = 0; i < 101; i++) {

        LOG(INFO) << i;

        printer.print(i, mesh.getCells(), mesh.getMyCells(), gas, size, rank);

        for (int j = 0; j < rep; ++j) {
            LOG(INFO) << "transfer";
            transfer->move(mesh.getFacets(), mesh.getCells(), mesh.getMyCells(), gas);
        }

        LOG(INFO) << "integral";
        integral.collide(rep*mesh.getTimeStep(), mesh.getCells(), mesh.getMyCells(), gas);

        for (int j = 0; j < rep; ++j) {
            LOG(INFO) << "transfer";
            transfer->move(mesh.getFacets(), mesh.getCells(), mesh.getMyCells(), gas);
        }
        
    }

    MPI_Finalize();

    return 0;
}
