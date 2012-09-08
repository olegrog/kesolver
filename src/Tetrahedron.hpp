#ifndef TETRAHEDRON_H_
#define TETRAHEDRON_H_

#include "Polygon.hpp"

class Tetrahedron : public Polygon {
	public:
 		void calculateVolume();			// вычисление объема тетраэдра
		void calculateCenter();			// вычисление центра тетраэдра

		Tetrahedron(); 		// конструктор класса
};

inline V3d tetrahedronCenter(V3d p1, V3d p2, V3d p3, V3d p4) {
	return (p1 + p2 + p3 + p4) / 4.;
}

inline double tetrahedronVolume(V3d p1, V3d p2, V3d p3, V3d p4) {
	return std::abs(det(p2 - p1, p3 - p1, p4 - p1)) / 6.;
}

#endif /*TETRAHIDRON_H_*/

