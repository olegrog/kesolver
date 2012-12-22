#include <iostream>
#include "Transfer1.hpp"
#include "Constructors.hpp"

void Transfer1::move(MeshMpi& mesh, const Gas& gas)
{
	mesh.newStep();

    for (size_t i = 0; i < mesh.getFlowingFacets().size(); i++) {
        mesh.getFlowingFacets()[i]->transfer(mesh.getAllCells(), gas);
    }
    for (size_t i = 0; i < mesh.getMyCells().size(); i++) 
        mesh.getMyCells()[i]->f().equategf();
}

