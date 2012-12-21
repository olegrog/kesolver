#ifndef POLYGON_H_
#define POLYGON_H_

#include <string>
#include <vector>

#include "SpeedFunction.hpp"
#include "v.hpp"

class Polygon {
    protected:

        int numberOfVertex;
        int numberOfEdges;
        int numberOfFacets;

        double V;
        double Lmin;
        V3d center;             // necessary just for second order
        T3d dd;

        std::vector<V3d> vertex;
        std::vector<int> neigbors;

        SpeedFunction function;

    public:

        int rank;
        std::string phys_name;

        void calculateLength();
        virtual void calculateVolume() = 0;
        virtual void calculateCenter() = 0;

        void setVertexes(const std::vector<V3d>& vertexes) { vertex = vertexes; }
        void setNeigbors(const std::vector<int>& neigbors_) { neigbors = neigbors_; }

        const std::vector<int>& getNeigbors() const { return neigbors; }

        double getLMin() const { return Lmin; }
        int getNumberOfVertexes() const { return  numberOfVertex; }
        
        void setRank(int rank_) { rank = rank_; }
        int getRank() const { return rank; }

        void setPhysicalName(const std::string& name) { phys_name = name; }
        const std::string& getPhysicalName() const { return phys_name; }

        double getVolume() const { return V; }
        V3d getCenter() const { return center; }
        T3d& getDD() { return dd; } 

        void inverseDD();

        void prepareForNextStep();
        void findGradient();

        SpeedFunction& f() { return function; }

};

#endif /* POLYGON_H_ */

