#ifndef _MESH_HPP_
#define _MESH_HPP_

#include <vector>

#include "mesh/MeshBase.hpp"

class Mesh {
    public:
        typedef typename MeshBase::Cells  Cells;
        typedef typename MeshBase::Facets Facets;

        virtual Cells&  getAllCells() = 0; 
        virtual Cells&  getFlowingCells() = 0;
        virtual Cells&  getMyCells() = 0;     

        virtual Facets& getAllFacets() = 0; 
        virtual Facets& getFlowingFacets() = 0;

        typedef std::vector<int> Ints;

        virtual const Ints& getMyCellIndexes() const = 0;

        virtual double getTimeStep() const = 0;

        virtual void newStep() = 0;

        virtual int getRank() const = 0;
        virtual int getSize() const = 0;

        virtual void barrier() const = 0;

};

#endif // _MESH_HPP_

