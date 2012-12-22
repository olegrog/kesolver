#ifndef TRANSFER_H
#define TRANSFER_H

#include <vector>

#include "Gas.hpp"
#include "Mesh.hpp"

class Transfer {
    public:
        virtual void move(Mesh& mesh, const Gas& gas) = 0;
};

#endif /*TRANSFER_H*/
