#include <iostream>
#include "Transfer1.hpp"
#include "Constructors.hpp"

void Transfer1::move(const MeshMpi& mesh, const Gas& gas)
{
	data_exchanger.swap();

    for(size_t i = 0; i < facets.size(); ++i)
        facets[i]->transfer(spacemesh, gas);
    for(size_t i = 0; i < mypolys.size(); ++i) 
        spacemesh[mypolys[i]]->f().equategf();

}

