#ifndef _TRANSFER2_HPP_
#define _TRANSFER2_HPP_

#include "Transfer.hpp"

class Transfer2 : public Transfer {
	public:
		Transfer2(Mesh& mesh);
		void move(Mesh& mesh,const Gas& gas);
};

#endif /*_TRANSFER2_HPP_*/
