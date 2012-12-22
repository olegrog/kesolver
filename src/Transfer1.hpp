#ifndef TRANSFER1_HPP_
#define TRANSFER1_HPP_

#include "Transfer.hpp"

class Transfer1 : public Transfer {
    public:
        void move(MeshMpi& mesh, const Gas& gas);
};

#endif /*TRANSFER1_HPP_*/
