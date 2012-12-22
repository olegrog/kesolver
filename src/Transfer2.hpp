#ifndef _TRANSFER2_HPP_
#define _TRANSFER2_HPP_

#include "Transfer.hpp"

class Transfer2 : public Transfer {
	public:
		Transfer2(MeshMpi& mesh);
		void move(MeshMpi& mesh,const Gas& gas);
};

#endif /*_TRANSFER2_HPP_*/
