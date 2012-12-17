#include <iostream>
#include <vector>
#include <string>

#include "mpi.h"

#include "property_tree/property_tree.hpp"

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

    LABEL
    std::cout << "gas.size() = " << gas_p->size() << std::endl;

    std::vector<Polygon*>        spacemesh;                   
    std::vector<PhysicalFacet*>  facets;  
    ElementsConstructor(prop_tree["mesh"], spacemesh, facets, gas);

    double curnt = prop_tree["curnt_limit"].asDouble();
    std::cout << "curnt = " << curnt << std::endl;

    double time_step = findTimeStep(spacemesh, gas, curnt);
    std::cout << "time_step = " << time_step << std::endl;
    
    for (std::vector<PhysicalFacet*>::iterator pp = facets.begin();
            pp != facets.end(); ++pp)
        (*pp)->findMultInOut(time_step, spacemesh);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::vector<int> mypolys;  
    MypolysConstructor(rank, spacemesh, mypolys);

    std::cout << "rank = " << rank << " facets.size() = "   << facets.size() << 
        " mypolys.size() = " << mypolys.size() << std::endl;

    Transfer* transfer;
    int order = prop_tree.isMember("order") ? prop_tree["order"].asInt() : 1;
    std::cout << "order = " << order << std::endl;

    if (order == 2)
        transfer = new Transfer2(facets, spacemesh, mypolys);
    else
        transfer = new Transfer1();
    transfer->init(prop_tree, gas, facets, spacemesh, mypolys, rank);

    Integral integral = IntegralConstructor(prop_tree);

    int rep = prop_tree.isMember("rep") ? prop_tree["rep"].asInt() : 1;
    std::cout << "rep = " << rep << std::endl;

    Printer printer(prop_tree);
    for (int i = 0; i < 2000000; i++) {

        std::cout << i << std::endl;

        printer.print(i, spacemesh, mypolys, gas, size, rank);

        for (int j = 0; j < rep; ++j) {
            std::cout << "transfer" << std::endl;
            transfer->move(facets, spacemesh, mypolys, gas);
        }

        std::cout << "integral" << std::endl;
        integral.collide(rep*time_step, spacemesh, mypolys, gas);

        for (int j = 0; j < rep; ++j) {
            std::cout << "transfer" << std::endl;
            transfer->move(facets, spacemesh, mypolys, gas);
        }
        
    }

    MPI_Finalize();

    return 0;
}
