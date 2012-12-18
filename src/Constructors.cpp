#include <algorithm>
#include <stdexcept>
#include <sstream>

#include "base64/base64.hpp"

#include "Constructors.hpp"

#include "auxiliary.hpp"

#include "Tetrahedron.hpp"
#include "Prism.hpp"
#include "Hexahedron.hpp"

#include "GateFacet.hpp"
#include "MirrorFacet.hpp"
#include "WallMaxwellFacet.hpp"
#include "MaxwellFacet.hpp"

#include "ximesh.hpp"
#include "ximesh_mixture.hpp"
#include "ximesh_mix.hpp"
#include "ximesh_rot.hpp"

#include "ci_simple.hpp"
#include "ci_mixture.hpp"
#include "ci_mix.hpp"

const std::vector<double> readMasses(const PropertyTree& tree) {
    std::vector<double> masses;
    if (tree.isMember("mass_ratio")) {
        const double m1   = 1;
        const double m2   = tree["mass_ratio"].asDouble();
        const double m0   = 0.5 * (m1 + m2);
        masses.push_back(m1 / m0);
        masses.push_back(m2 / m0);
    }
    else {
        std::istringstream ss(tree["masses"].asString());
        double mass;
        while (ss >> mass)
            masses.push_back(mass);
    }
    return masses;
}

template <Symmetry symmetry, Volume volume, TimeScheme scheme>
Gas* gasMixScheme(const PropertyTree& tree, 
                  const int rad, const double cut,
                  const std::vector<double>& masses)
{
    typedef typename SymmetryTrait<symmetry>::Vm Vm;
    Vm v;
    if (tree.isMember("v")) 
        v = strTo<Vm>(tree["v"].asString());
    else 
        v = 0.;
    std::vector<double> adds;
    if (tree.isMember("adds")) {
        std::istringstream ss(tree["adds"].asString());
        double add;
        while (ss >> add)
            adds.push_back(add);
    }
    else {
        for (size_t i = 0; i < masses.size(); ++i)
            adds.push_back(0.0);
    }
    typedef GasTemplate<symmetry, XiMeshMix, ColliderMix<symmetry, volume, scheme> > GasType;
    return new GasType(XiMeshMix<symmetry>(rad, cut, masses, v, adds));
}

template <Symmetry symmetry, Volume volume>
Gas* gasMixVolume(const PropertyTree& tree, 
                  const int rad, const double cut,
                  const std::vector<double>& masses)
{
    std::string scheme;
    if (tree.isMember("time_scheme")) 
        scheme = tree["time_scheme"].asString();
    else
        scheme = "Continues";

    if      (scheme == "Euler") 
        return gasMixScheme<symmetry, volume, Euler>     (tree, rad, cut, masses);
    else if (scheme == "Continues")
        return gasMixScheme<symmetry, volume, Continues> (tree, rad, cut, masses);
    else { return NULL; /* TODO raise */ }
}

template <Symmetry symmetry>
Gas* gasSymmetry(const PropertyTree& tree, const std::string& type, 
                 const int rad, const double cut)
{
    LABEL
    if (type == "Simple") {
        LABEL
        typedef typename SymmetryTrait<symmetry>::Vm Vm;
        Vm v;
        if (tree.isMember("v")) 
            v = strTo<Vm>(tree["v"].asString());
        else 
            v = 0.;
        std::cout << "v = " << v << std::endl;
        typedef GasTemplate<symmetry, XiMesh, ColliderSimple<symmetry> > GasType;
        return new GasType(XiMesh<symmetry>(rad, cut, v));
    }
    else if (type == "Mixture") {
        typedef GasTemplate<symmetry, XiMeshMixture, ColliderMixture<symmetry> > GasType;
        const std::vector<double> masses = readMasses(tree);

        const double a = cut / rad * std::sqrt(masses[0]);
        std::vector<int> rads;
        for (size_t i = 0; i < masses.size(); ++i)
            rads.push_back(static_cast<int>(rad * std::sqrt(masses[i] / masses[0])));
        return new GasType(XiMeshMixture<symmetry>(a, rads, masses));
    }
    else if (type == "Mix") {
        const std::vector<double> masses = readMasses(tree);
        const std::string volume = tree["volume"].asString();
        if      (volume == "Wide") 
            return gasMixVolume<symmetry, Wide>      (tree, rad, cut, masses);
        else if (volume == "Symmetric")
            return gasMixVolume<symmetry, Symmetric> (tree, rad, cut, masses);
        else if (volume == "Tight")
            return gasMixVolume<symmetry, Tight>     (tree, rad, cut, masses);
        else if (volume == "Tight2")
            return gasMixVolume<symmetry, Tight2>    (tree, rad, cut, masses);
        else if (volume == "Grad")
            return gasMixVolume<symmetry, Grad>      (tree, rad, cut, masses);
        else { return NULL; /* TODO raise */ }
    }
    else { return NULL; /* TODO raise */ }
}

void GasConstructor(const PropertyTree& tree, Gas** gas_pp)
{
    const std::string type = tree["type"].asString();
    std::cout << "type = " << type << std::endl;

    int    rad = tree["rad"].asInt();
    double cut = tree["cut"].asDouble();
    std::cout << "rad, cut = " << rad << ' ' << cut << std::endl;

    std::string symm;
    if (tree.isMember("symmetry"))
        symm = tree["symmetry"].asString();
    else 
        symm = "Cartesian";
    std::cout << "symm = " << symm << std::endl;

    Gas* gas_p;
    if (symm == "Cylindrical")
        gas_p = gasSymmetry<Cylindrical>(tree, type, rad, cut);
    else if (symm == "Cartesian")
        gas_p = gasSymmetry<Cartesian>(tree, type, rad, cut);
    else { gas_p = NULL; /* TODO raise */ }

    *gas_pp = gas_p;
}


Polygon* createPolygon(const PropertyTree& celldata) {
    int type = celldata["type"].asInt();
    if      (type == 4) return new Tetrahedron();
    else if (type == 5) return new Hexahedron();
    else if (type == 6) return new Prism();
    else return 0;
}

PhysicalFacet* createFacet(const std::string& type,
                           const PropertyTree& physdata, 
                           const Gas& gas)
{
    if (type == "diffusion")  {
        return new WallMaxwellFacet(physdata, gas);
    }
    else if (type == "maxwell")  {
        return new MaxwellFacet(physdata, gas);
    }
    else if (type == "mirror") {
        return new MirrorFacet((Axis)physdata["axis"].asInt());
    }
    else if (type == "gate") {
        return new GateFacet();
    }
    else {
        // TODO: raise exception
        return 0;
    }
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
    std::string phys_name = celldata["phys_name"].asString();

    polygon->setRank(rank);
    polygon->setPhysicalName(phys_name);

    return polygon;
}

PhysicalFacet* constructFacet(const PropertyTree& facetdata, 
                              const PropertyTree& bcsdata, 
                              const std::vector<V3d>& nodes, 
                              const Gas& gas)
{
    std::string phys_name = facetdata["phys_name"].asString();
    std::cout << "facet: phys_name = " << phys_name << std::endl;
    PhysicalFacet* facet = bcsdata.isMember(phys_name) ? 
                           createFacet(bcsdata[phys_name]["type"].asString(),
                                       bcsdata[phys_name], gas) :
                           createFacet("gate", facetdata, gas);
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
    std::cout << "nodes_num = " << nodes_num << std::endl;
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
    std::cout << "cells.size() = " << polygons.size() << std::endl;

    // construct facets
    const PropertyTree& facetsdata = meshdata["facets"];
    for (size_t i = 0; i < facetsdata.size(); ++i)
        facets.push_back(constructFacet(facetsdata[i], bcsdata, nodes, gas));
    std::sort(facets.begin(), facets.end(), lessFacet);
    std::cout << "facets.size() = " << facets.size() << std::endl;

}

void MypolysConstructor(int rank,
                        const std::vector<Polygon*>& polygons,
                        std::vector<int>& mypolys)
{
    for (size_t i = 0; i < polygons.size(); ++i) 
        if (polygons[i]->getRank() == rank) {
            mypolys.push_back(i);
        }
}

double findTimeStep(const std::vector<Polygon*>& spacemesh,
                    const Gas& gas,
                    double curnt)
{
    double h = std::numeric_limits<double>::max();
    for (size_t i = 0; i < spacemesh.size(); ++i) 
        if (spacemesh[i]->getLMin() < h) 
            h = spacemesh[i]->getLMin();
    std::cout << "Lmin = " << h << std::endl;

    double xi_max = 0.0;
    for (size_t i = 0; i < gas.size(); ++i) {
        xi_max = std::max(xi_max, gas.dot(i, V3d(1., 0., 0.)));
    }
    std::cout << "xi_max = " << xi_max << ' ' << gas.cut() << std::endl;

    return curnt * h / gas.cut();
}


void GivePolygonMemoryAndInit(const PropertyTree& tree, const Gas& gas, Polygon* polygon)
{
    const std::string name = polygon->getPhysicalName();
    std::cout << "give and init: " << name << std::endl;
    const PropertyTree& data = tree["initial_conditions"][name];
    polygon->f().f(gas.maxwell(data));
}

Integral IntegralConstructor(const PropertyTree& tree)
{
    int intorder;
    if (tree.isMember("order"))
        intorder = tree["order"].asInt();
    else
        intorder = 1;
    std::cout << "intorder = " << intorder << std::endl;

    bool is_free_molecular = !(tree["enable"].asBool());

    int power = tree["power"].asInt();

    SimpleSection* section;
    if (tree.isMember("section"))
    {
        std::string str = tree["section"].asString();
        if (str == "Simple") {
            throw std::invalid_argument("Simple section");
        }
        else if (str == "HS") {
            std::vector<double> ds;
            std::istringstream ss(tree["ds"].asString());
            double d;
            while (ss >> d)
                ds.push_back(d);
            section = new HSSection(ds);
            std::cout << "HSSection" << std::endl;

        }
        else if (str == "LJ") {
            std::vector<double> ds, es;
            std::istringstream ss(tree["ds"].asString());
            double x;
            while (ss >> x)
                ds.push_back(x);
            std::istringstream ss2(tree["es"].asString());
            while (ss2 >> x)
                es.push_back(x);
            std::string filename = tree["file"].asString();
            section = new LJSection(ds, es, filename);
            std::cout << "LJSection" << std::endl;
        }
        else if (str == "Ma") {
            std::vector<double> as;
            std::istringstream ss(tree["ds"].asString());
            double x;
            while (ss >> x)
                as.push_back(x);
            std::string filename = tree["file"].asString();
            section = new MaSection(as, filename);
            std::cout << "MaSection" << std::endl;
        }
        else
            throw std::invalid_argument("Unknown section");
    }
    else {
        section = new SimpleSection;
        std::cout << "SimpleSection" << std::endl;
    }

    return Integral(power, intorder, is_free_molecular, section);
}

