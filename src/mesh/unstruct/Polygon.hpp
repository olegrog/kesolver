#ifndef POLYGON_H_
#define POLYGON_H_

#include <string>
#include <vector>

#include "mesh/CellBase.hpp"

#include "Gas.hpp"
#include "SpeedFunction.hpp"
#include "v.hpp"

class Polygon : public CellBase {
    protected:

        int numberOfVertex;
        int numberOfEdges;
        int numberOfFacets;

        double V;
        double Lmin;
        V3d center;             // necessary just for second order
        T3d dd;

        std::vector<V3d> vertex;

        SpeedFunction function;

    public:
        void calculateLength();
        virtual void calculateVolume() = 0;
        virtual void calculateCenter() = 0;

        void setVertexes(const std::vector<V3d>& vertexes) { vertex = vertexes; }

        double getLMin() const { return Lmin; }
        int getNumberOfVertexes() const { return  numberOfVertex; }
        
        void setRank(int rank_) { rank = rank_; }
        int getRank() const { return rank; }

        double getVolume() const { return V; }
        V3d getCenter() const { return center; }
        T3d& getDD() { return dd; } 

        void inverseDD();

        void prepareForNextStep();
        void findGradientAndPhi(const Gas& gas);

        SpeedFunction& f() { return function; }

};

#endif /* POLYGON_H_ */

