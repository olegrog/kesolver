#ifndef _CELLBASE_H_
#define _CELLBASE_H_

#include <vector>

// TOREMOVE:
#include "SpeedFunction.hpp"

class CellBase {
    public:
        typedef std::vector<int> Ints;

        void setIndex(int i) { index = i; }
        int getIndex() const { return index; }

        void setNeigbors(const Ints& ns) { neigbors = ns; }
        const Ints& getNeigbors() const { return neigbors; }

        void setRank(int r) { rank = r; }
        int getRank() const { return rank; }

        void setPhysicalName(const std::string& name) { phys_name = name; }
        const std::string& getPhysicalName() const { return phys_name; }

        void setExcludedPoints(double points) { excluded_points = points; }
        double getExcludedPoints() const { return excluded_points; }

        CellBase() : excluded_points(0) {}

        // TOREMOVE:
        virtual SpeedFunction& f() = 0;

    protected:
        int index;
        Ints neigbors;
        std::string phys_name;
        int rank;
        double excluded_points;
};

#endif /* _CELLBASE_H_ */

