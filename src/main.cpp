#include <iostream>
#include <vector>
#include <string>
#include <limits>
#include <stdexcept>
#include "mpi.h"

#include "Polygon.hpp"
#include "Loader.hpp"
#include "Constructors.hpp"
#include "Printer.hpp"
#include "Integral.hpp"
#include "Transfer.hpp"
#include "Transfer1.hpp"
#include "Transfer2.hpp"
#include "DataExchanger.hpp"
#include "Tetrahedron.hpp"
#include "PhysicalFacet.hpp"
#include "Gas.hpp"

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    std::string kinfile_name(argv[1]);
    Loader loader(kinfile_name);                    

    Gas* gas_p;
    GasConstructor(loader, &gas_p);
    LABEL
    std::cout << "gas.size() = " << gas_p->size() << std::endl;
    Gas& gas = *gas_p;

    std::vector<Polygon*> spacemesh;                   
    PolygonsConstructor(loader, spacemesh);

    double curnt = loader.getData<double>("curnt", "value");
    std::cout << "curnt = " << curnt << std::endl;
    double time_step = findTimeStep(spacemesh, gas, curnt);
    std::cout << "time_step = " << time_step << std::endl;

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::vector<PhysicalFacet*> facets;  
    FacetsConstructor(loader, spacemesh, facets, gas, time_step, rank);

    std::vector<int> mypolys;  
    MypolysConstructor(rank, spacemesh, mypolys);

    std::cout << "rank = " << rank << " facets.size() = "   << facets.size() << 
        " mypolys.size() = " << mypolys.size() << std::endl;

    Transfer* transfer;
    int order;
    try {
        order = loader.getData<int>("transfer", "order");
    }
    catch (std::invalid_argument) {
        std::cout << "catch" << std::endl;
        order = 1;
    }
    std::cout << "order = " << order << std::endl;
    if (order == 2)
        transfer = new Transfer2(facets, spacemesh, mypolys);
    else
        transfer = new Transfer1();
    transfer->init(loader, gas, facets, spacemesh, mypolys, rank);

    Integral integral = IntegralConstructor(loader);

    int rep;
    try {
        rep = loader.getData<int>("transfer", "rep");
    }
    catch (std::invalid_argument) {
        std::cout << "catch" << std::endl;
        rep = 1;
    }
    std::cout << "rep = " << rep << std::endl;

    Printer printer(loader);
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
