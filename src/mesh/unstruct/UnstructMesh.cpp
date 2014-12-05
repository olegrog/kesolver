#include <limits>

#include "UnstructMesh.hpp"

#include "Tetrahedron.hpp"
#include "Prism.hpp"
#include "Hexahedron.hpp"

#include "FacetFactory.hpp"

#include "base64/base64.hpp"
#include "logger/logger.hpp"

Polygon* createPolygon(const PropertyTree& celldata) {
    int type = celldata["type"].asInt();
    if      (type == 4) return new Tetrahedron();
    else if (type == 5) return new Hexahedron();
    else if (type == 6) return new Prism();
    else return 0;
}

template <typename Element>
void constructElement(const PropertyTree& elemdata, 
                      const std::vector<V3d>& nodes, 
                      Element* elem)
{
    const PropertyTree& ns = elemdata["nodes"];
    std::vector<V3d> vertexes(ns.size());
    for (size_t i = 0; i < vertexes.size(); ++i) 
        vertexes[i] = nodes[ns[i].asInt()];

    const PropertyTree& nb = elemdata["neigbors"];
    std::vector<int> neigbors(nb.size());
    for (size_t i = 0; i < neigbors.size(); ++i) 
        neigbors[i] = nb[i].asInt();

    elem->setVertexes(vertexes);
    elem->setNeigbors(neigbors);
}

Polygon* constructPolygon(const PropertyTree& celldata, 
                      const std::vector<V3d>& nodes)
{
    Polygon* polygon = createPolygon(celldata);

    constructElement(celldata, nodes, polygon);

    polygon->calculateLength();
    polygon->calculateVolume();
    polygon->calculateCenter();

    int rank = celldata["part_index"].asInt();
    int index = celldata["ord_index"].asInt();
    std::string phys_name = celldata["phys_name"].asString();

    polygon->setRank(rank);
    polygon->setIndex(index);
    polygon->setPhysicalName(phys_name);

    return polygon;
}

PhysicalFacet* constructFacet(const PropertyTree& facetdata, 
                              const PropertyTree& bcsdata, 
                              const std::vector<V3d>& nodes, 
                              const Gas& gas)
{
    std::string phys_name = facetdata["phys_name"].asString();
    std::string  facet_name = "gate";
    PropertyTree facet_cond = facetdata;
    if (bcsdata.isMember(phys_name)) {
        facet_name = bcsdata[phys_name]["type"].asString();
        facet_cond = bcsdata[phys_name];
    }
    PhysicalFacet* facet = FacetFactory::instance().create(facet_name);
    facet->init(facet_cond, gas);
    constructElement(facetdata, nodes, facet);
    facet->setType(facetdata["type"].asInt());
    facet->findNormalAndSquare();

    return facet;
}


bool lessFacet(PhysicalFacet* f1, PhysicalFacet* f2) {
    return f1->order() < f2->order();
}

void ElementsConstructor(const PropertyTree& tree,
                         std::vector<Polygon*>& polygons,
                         std::vector<PhysicalFacet*>& facets, 
                         const Gas& gas)
{
    const PropertyTree& meshdata = tree["mesh"];
    const PropertyTree& bcsdata  = tree["boundary_conditions"];

    // construct nodes
    std::vector<unsigned char> nodesbytes = base64::decode(meshdata["nodes"].asString());
    const double* nodesdoubles = reinterpret_cast<const double*>(&nodesbytes.front());
    size_t nodes_num = meshdata["nodes_num"].asInt();
    LOG(INFO) << "nodes_num = " << nodes_num;
    std::vector<V3d> nodes(nodes_num);
    for (size_t i = 0, j = 0; i < nodes.size(); ++i) {
        double x = nodesdoubles[j++];
        double y = nodesdoubles[j++];
        double z = nodesdoubles[j++];
//      std::cout << x << ' ' << y << ' ' << z << std::endl;
        nodes[i] = V3d(x, y, z);
    }

    // construct cells
    const PropertyTree& cellsdata = meshdata["cells"];
    for (size_t i = 0; i < cellsdata.size(); ++i) 
        polygons.push_back(constructPolygon(cellsdata[i], nodes));
    LOG(INFO) << "cells.size() = " << polygons.size();

    // construct facets
    const PropertyTree& facetsdata = meshdata["facets"];
    for (size_t i = 0; i < facetsdata.size(); ++i)
        facets.push_back(constructFacet(facetsdata[i], bcsdata, nodes, gas));
    std::sort(facets.begin(), facets.end(), lessFacet);
    LOG(INFO) << "facets.size() = " << facets.size();

}

double findTimeStep(const std::vector<Polygon*>& spacemesh,
                    const Gas& gas,
                    double curnt)
{
    double h = std::numeric_limits<double>::max();
    for (size_t i = 0; i < spacemesh.size(); ++i) 
        if (spacemesh[i]->getLMin() < h) 
            h = spacemesh[i]->getLMin();
    LOG(INFO) << "Lmin = " << h;

    double xi_max = 0.0;
    for (size_t i = 0; i < gas.size(); ++i) {
        xi_max = std::max(xi_max, gas.dot(i, V3d(1., 0., 0.)));
    }
    LOG(INFO) << "xi_max = " << xi_max << ' ' << gas.cut();

    return curnt * h / gas.cut();
}

UnstructMesh::UnstructMesh(const PropertyTree& tree,
                           const Gas& gas)
{
    ElementsConstructor(tree, cells, facets, gas);

    double curnt = tree["curnt_limit"].asDouble();
    LOG(INFO) << "curnt = " << curnt;

    time_step = findTimeStep(cells, gas, curnt);
    LOG(INFO) << "time_step = " << time_step;
    
    for (std::vector<PhysicalFacet*>::iterator pp = facets.begin();
            pp != facets.end(); ++pp)
        (*pp)->findMultInOut(time_step, cells);
}

