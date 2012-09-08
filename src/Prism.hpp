#ifndef PRISM_H_
#define PRISM_H_

#include "Polygon.hpp"

class Prism : public Polygon {
	public:
		void calculateVolume();			// вычисление объема
		void calculateCenter();			// вычисление центра

		Prism(); 		// конструктор класса
};

#endif /* PRISM_H_ */
