#include <iostream>
#include <vector>
#include <string>

#include "mpi.h"

#include "property_tree/property_tree.hpp"
#include "logger/logger.hpp"

#include "Constructors.hpp"
#include "Polygon.hpp"
#include "PhysicalFacet.hpp"

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

    std::vector<Polygon*>        spacemesh;                   
    std::vector<PhysicalFacet*>  facets;  
    ElementsConstructor(prop_tree, spacemesh, facets, gas);

    double curnt = prop_tree["curnt_limit"].asDouble();
    LOG(INFO) << "curnt = " << curnt;

    double time_step = findTimeStep(spacemesh, gas, curnt);
    LOG(INFO) << "time_step = " << time_step;
    
    for (std::vector<PhysicalFacet*>::iterator pp = facets.begin();
            pp != facets.end(); ++pp)
        (*pp)->findMultInOut(time_step, spacemesh);

    LABEL

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::vector<int> mypolys;  
    MypolysConstructor(rank, spacemesh, mypolys);

    LOG(INFO) << "rank = " << rank << " facets.size() = " << facets.size() << 
        " mypolys.size() = " << mypolys.size();

    Transfer* transfer;
    int order = prop_tree.isMember("order") ? prop_tree["order"].asInt() : 1;
    LOG(INFO) << "order = " << order;

    if (order == 2)
        transfer = new Transfer2(facets, spacemesh, mypolys);
    else
        transfer = new Transfer1();
    transfer->init(prop_tree, gas, facets, spacemesh, mypolys, rank);

    Integral integral = IntegralConstructor(prop_tree["integral"]);

    int rep = prop_tree["transfer"].isMember("rep") ? prop_tree["transfer"]["rep"].asInt() : 1;
    LOG(INFO) << "rep = " << rep;

    Printer printer(prop_tree["printer"]);
    for (int i = 0; i < 101; i++) {

        LOG(INFO) << i;

        printer.print(i, spacemesh, mypolys, gas, size, rank);

        for (int j = 0; j < rep; ++j) {
            LOG(INFO) << "transfer";
            transfer->move(facets, spacemesh, mypolys, gas);
        }

        LOG(INFO) << "integral";
        integral.collide(rep*time_step, spacemesh, mypolys, gas);

        for (int j = 0; j < rep; ++j) {
            LOG(INFO) << "transfer";
            transfer->move(facets, spacemesh, mypolys, gas);
        }
        
    }

    MPI_Finalize();

    return 0;
}
