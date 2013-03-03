#include <iostream>
#include <vector>
#include <string>

#include "property_tree/property_tree.hpp"
#include "logger/logger.hpp"

#include "Constructors.hpp"

#include "mesh/unstruct/UnstructMesh.hpp"
#include "Mesh.hpp"
#include "MeshSingle.hpp"

#ifdef MPI_ON
#include "MeshMpi.hpp"
#endif

#include "Transfer.hpp"
#include "Transfer1.hpp"
#include "Transfer2.hpp"

#include "Integral.hpp"

#include "Printer.hpp"

int main(int argc, char** argv)
{
    PropertyTree prop_tree(argv[1]);

    Gas* gas_p;
    GasConstructor(prop_tree["gas"], &gas_p);
    Gas& gas = *gas_p;

    LOG(INFO) << "gas.size() = " << gas_p->size();

    MeshBase* mesh_base_p = new UnstructMesh(prop_tree, gas);
#ifdef MPI_ON
    MeshMpi*     mesh_p   = new MeshMpi(argc, argv, mesh_base_p);
#else
    MeshSingle*  mesh_p   = new MeshSingle(mesh_base_p);
#endif
    Mesh& mesh = *mesh_p;

    LOG(INFO) << "facets.size() = "    << mesh.getAllFacets().size() 
              << " ffacets.size() = "  << mesh.getFlowingFacets().size();

    LOG(INFO) << "cells.size() = "     << mesh.getAllCells().size() 
              << " fcells.size() = "   << mesh.getFlowingCells().size()
              << " mycells.size() = "  << mesh.getMyCells().size();

	for (size_t i = 0; i < mesh.getFlowingCells().size(); ++i) {
		GivePolygonMemoryAndInit(prop_tree, gas, mesh.getFlowingCells()[i]);
    }

    Transfer* transfer;
    int order = prop_tree["transfer"].isMember("order") ? 
                prop_tree["transfer"]["order"].asInt() : 1;
    LOG(INFO) << "order = " << order;

    if (order == 2)
    {
	    for (size_t i = 0; i < mesh.getFlowingCells().size(); ++i) {
		    mesh.getFlowingCells()[i]->f().giveMemoryToGradient();
        }
        transfer = new Transfer2(mesh);
    }
    else
        transfer = new Transfer1();

    Integral integral = IntegralConstructor(prop_tree["integral"]);

    int rep = prop_tree["transfer"].isMember("rep") ? 
              prop_tree["transfer"]["rep"].asInt() : 1;
    LOG(INFO) << "rep = " << rep;

    Printer printer(prop_tree["printer"]);
    for (int i = 0; i < 200001; i++) {

        LOG(INFO) << i;

        printer.print(i, mesh, gas);

        for (int j = 0; j < rep; ++j) {
            LOG(INFO) << "transfer";
            transfer->move(mesh, gas);
        }

        LOG(INFO) << "integral";
        double ts = rep*mesh.getTimeStep();
        LABEL
        integral.collide(ts, mesh.getMyCells(), gas);

        for (int j = 0; j < rep; ++j) {
            LOG(INFO) << "transfer";
            transfer->move(mesh, gas);
        }
        
    }

    return 0;
}
