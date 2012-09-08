#ifndef HEXAHEDRON_H_
#define HEXAHEDRON_H_

#include "Polygon.hpp"

class Hexahedron : public Polygon {
	public:
		void calculateVolume();			// вычисление объема
		void calculateCenter();			// вычисление центра

		Hexahedron(); 		// конструктор класса
};

#endif /* HEXAHEDRON_H_ */
