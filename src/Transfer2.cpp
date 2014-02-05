#include <iostream>
#include "Transfer2.hpp"

Transfer2::Transfer2(Mesh& mesh)
{
	for (size_t i = 0; i < mesh.getFlowingFacets().size(); i++) 
		mesh.getFlowingFacets()[i]->calculateDistance(mesh.getAllCells());
	for (size_t i = 0; i < mesh.getFlowingCells().size(); i++) 
		mesh.getFlowingCells()[i]->inverseDD();
}

void Transfer2::move(Mesh& mesh, const Gas& gas)
{
	mesh.newStep();

    for (size_t i = 0; i < mesh.getFlowingCells().size(); i++) {
        mesh.getFlowingCells()[i]->prepareForNextStep();
    }
    for (size_t i = 0; i < mesh.getFlowingFacets().size(); i++) {
        mesh.getFlowingFacets()[i]->findGradient(mesh.getAllCells(), gas);
    }
    for (size_t i = 0; i < mesh.getFlowingCells().size(); i++) {
        mesh.getFlowingCells()[i]->findGradientAndPhi(gas);
    }
    for (size_t i = 0; i < mesh.getFlowingFacets().size(); i++) {
        mesh.getFlowingFacets()[i]->findPhi(mesh.getAllCells(), gas);
    }
    for (size_t i = 0; i < mesh.getFlowingCells().size(); i++) {
        mesh.getFlowingCells()[i]->f().equatefg();
    }
    for (size_t i = 0; i < mesh.getMyFacets().size(); i++) {
        mesh.getMyFacets()[i]->transfer2(mesh.getAllCells(), gas);
    }
    for (size_t i = 0; i < mesh.getMyCells().size(); i++) 
        mesh.getMyCells()[i]->f().equategf();

}

