#ifndef POLYGON_H_
#define POLYGON_H_

#include "SpeedFunction.hpp"
#include "v.hpp"

class Polygon {
	protected:

		int numberOfVertex;
		int numberOfEdges;
		int numberOfFacets;

		double V;
		double Lmin;
		V3d center;				// necessary just for second order
		T3d dd;

		std::vector<V3d> vertex;
		std::vector<int> neigbors;

		SpeedFunction function;

	public:

		int rank;
		int physical_index;

		void calculateLength();
		virtual void calculateVolume() = 0;
		virtual void calculateCenter() = 0;

		void setVertexCoordinates(const std::vector<V3d>& vertex_) { vertex = vertex_; }
		void setNeigbors(const std::vector<int>& neigbors_) { neigbors = neigbors_; }

		const std::vector<int>& getNeigbors() const { return neigbors; }

		double getLMin() const { return Lmin; }
		int getNumberOfVertexes() const { return  numberOfVertex; }
		
		void setRank(int rank_) { rank = rank_; }
		int getRank() const { return rank; }

		void setPhysicalIndex(int index) { physical_index = index; }
		int getPhysicalIndex() const { return physical_index; }

		double getVolume() const { return V; }
		V3d getCenter() const { return center; }
		T3d& getDD() { return dd; } 

		void inverseDD();

		void prepareForNextStep();
		void findGradient();

		SpeedFunction& f() { return function; }

		bool isActive() const { return function.size() > 0; }

};

#endif /* POLYGON_H_ */

