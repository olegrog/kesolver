#ifndef _SOLVER_HPP_
#define _SOLVER_HPP_

#include <limits>

#include "property_tree/property_tree.hpp"
#include "logger/logger.hpp"

#include "Constructors.hpp"

#include "mesh/unstruct/UnstructMesh.hpp"
#include "Mesh.hpp"

#include "Transfer.hpp"
#include "Transfer1.hpp"
#include "Transfer2.hpp"

#include "Integral.hpp"

#include "Printer.hpp"

template <typename MeshType>
class Solver {
    public:
        Solver(int argc, char** argv) {

            LABEL
            PropertyTree prop_tree(argv[1]);
            LABEL

            Gas* gas_p;
            GasConstructor(prop_tree["gas"], &gas_p);
            Gas& gas = *gas_p;

            LOG(INFO) << "gas.size() = " << gas_p->size();

            MeshBase* mesh_base_p = new UnstructMesh(prop_tree, gas);
            Mesh*     mesh_p   = new MeshType(argc, argv, mesh_base_p);
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

            int num_time_steps = prop_tree.isMember("num_steps") ?
                                 prop_tree["num_steps"].asInt() : std::numeric_limits<int>::max();

            for (int i = 0; i < num_time_steps; i++) {

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

        }
};


#endif
