#include <iostream>

#include "MeshMpi.hpp"
#include "solver.hpp"

int main(int argc, char** argv)
{
    Solver<MeshMpi>(argc, argv);

    return 0;
}
