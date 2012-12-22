#ifndef TRANSFER_H
#define TRANSFER_H

#include <vector>

#include "Gas.hpp"
#include "MeshMpi.hpp"

class Transfer {
    public:
        virtual void move(MeshMpi& mesh, const Gas& gas) = 0;
};

#endif /*TRANSFER_H*/
