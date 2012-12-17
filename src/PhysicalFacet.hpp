#ifndef PHYSICALFACET_H_
#define PHYSICALFACET_H_

#include <vector>

#include "Gas.hpp"
#include "Polygon.hpp"

class PhysicalFacet {
	protected:
		int type;

		double S;
		V3d n, center;

		double mult_in, mult_out;
        double dt;

		int numberOfVertex;
		std::vector<V3d> vertex;
		std::vector<int> polygon;
		V3d d_in, d_out;

		bool is_active;

	public:
		void setType(int type_);
		void setNeigbors(const std::vector<int>& neigbors_);
		void setVertexes(const std::vector<V3d>& vertex_);
		void findNormalAndSquare();
		void findMultInOut(double t, const std::vector<Polygon*>& spacemesh);
		double getSquare() const { return S; }
		V3d getCenter() const { return center; }

		void transfer(std::vector<Polygon*>& spacemesh, const Gas& gas);
		virtual void doTransfer(std::vector<Polygon*>& spacemesh,
				const Gas& gas) = 0;
		void transfer2(std::vector<Polygon*>& spacemesh, const Gas& gas);
		virtual void doTransfer2(std::vector<Polygon*>& spacemesh,
				const Gas& gas) = 0;

		void findPhi(const std::vector<Polygon*>& spacemesh, const Gas& gas);
		void doFindPhi(const std::vector<Polygon*>& spacemesh, const Gas& gas);

		virtual void calculateDistance(std::vector<Polygon*>& spacemesh) = 0;

		void findGradient(const std::vector<Polygon*>& spacemesh,
				const Gas& gas);
		virtual void doFindGradient(const std::vector<Polygon*>& spacemesh,
				const Gas& gas) = 0;

		static int getNumberOfVertexes(int type) {
			if (type == 2) return 3;
			else if (type == 3) return 4;
			else return 0;
		}
		virtual int order() const = 0;

		void activation(const std::vector<Polygon*>& spacemesh);

	private:

};

#endif /* PHYSICALFACET_H_ */
