#include <iostream>

#include "MeshSingle.hpp"
#include "solver.hpp"

int main(int argc, char** argv)
{
    Solver<MeshSingle>(argc, argv);

    return 0;
}
