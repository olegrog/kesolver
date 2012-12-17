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
    if      (symm == "Cylindrical")
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
    
    int rank           = elemdata["part_index"].asInt();
    int physical_index = elemdata["phys_index"].asInt();

    elem->setType(type);
    elem->setVertexes(vertexes);
    elem->setNeigbors(neigbors);

void constructPolygon(const PropertyTree& celldata, 
                      const std::vector<V3d>& nodes, 
                      Polygon* polygon)
{
    const PropertyTree& ns = celldata["nodes"];
    std::vector<V3d> vertexes(ns.size());
    for (size_t i = 0; i < vertexes.size(); ++i) 
        vertexes[i] = nodes[ns[i].asInt()];

    const PropertyTree& nb = celldata["nb"];
    std::vector<int> neigbors(nb.size());
    for (size_t i = 0; i < neigbors.size(); ++i) 
        neigbors[i] = nb[i].asInt();
    
    int rank           = celldata["part_index"].asInt();
    int physical_index = celldata["phys_index"].asInt();

    polygon->setVertexCoordinates(vertexes);
    polygon->calculateLength();
    polygon->calculateVolume();
    polygon->calculateCenter();
    polygon->setNeigbors(neigbors);
    polygon->setRank(rank);
    polygon->setPhysicalIndex(physical_index);
}

PhysicalFacet* createFacet(const PropertyTree& facetdata, const Gas& gas) {
    std::string type = facetdata["type"].asString();
    if (type == "diffusion")  {
        std::vector<std::string> newdata;
        std::copy(++data.begin(), data.end(), std::back_inserter(newdata));
        return new WallMaxwellFacet(gas.maxwell(newdata));
    }
    else if (type == "constant")  {
        std::vector<std::string> newdata;
        std::copy(++data.begin(), data.end(), std::back_inserter(newdata));
        return new MaxwellFacet(gas.maxwell(newdata));
    }
    else if (type == "mirror") {
        return new MirrorFacet((Axis)strTo<int>(data[1]));
    }
    else {
        return new GateFacet();
    }
}

bool lessFacet(PhysicalFacet* f1, PhysicalFacet* f2) {
    return f1->order() < f2->order();
}

void ElementsConstructor(const PropertyTree& tree,
                         std::vector<Polygon*>& polygons,
                         std::vector<PhysicalFacet*>& facets, 
                         const Gas& gas)
{
    const PropertyTree& cellsdata = tree["cells"];
    for (size_t i = 0; i < cellsdata.size(); ++i) 
        polygons.push_back(createPolygon(cellsdata[i]));

    std::vector<unsigned char> nodesbytes = base64::decode(tree["nodes"].asString());
    const double* nodesdoubles = reinterpret_cast<const double*>(&nodesbytes.front());
    std::vector<V3d> nodes(tree["nodes_num"].asInt());
    for (size_t i = 0, j = 0; i < nodes.size(); ++i) {
        double x = nodesdoubles[j++];
        double y = nodesdoubles[j++];
        double z = nodesdoubles[j++];
        nodes[i] = V3d(x, y, z);
    }

    for (size_t i = 0; i < cellsdata.size(); ++i) 
        constructPolygon(cellsdata[i], nodes, polygons, polygons[i]);

    const PropertyTree& facetsdata = tree["facets"];
    for (size_t i = 0; i < facetsdata.size(); ++i) {

        std::vector<int>::const_iterator p = (*pp).begin();
        int type = *p++;
        std::vector<V3d> vertexes(PhysicalFacet::getNumberOfVertexes(type));
        for (size_t i = 0; i < vertexes.size(); ++i) 
            vertexes[i] = loader.getNodes()[*p++];

        std::vector<int> neigbors(*p++);
        for (size_t i = 0; i < neigbors.size(); ++i) 
            neigbors[i] = *p++;

        int number_of_partitions = *p++;
        for (int i = 0; i < number_of_partitions; ++i) 
            p++;

        int key = *p;
        PhysicalFacet* facet = createFacet(loader.getPhysicalData(key), gas);
                
        facet->setType(type);
        facet->setVertex(vertexes);
        facet->setPolygonNumbers(neigbors);

        facets.push_back(facet);
    }
    std::sort(facets.begin(), facets.end(), lessFacet);

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

/*

void GivePolygonMemoryAndInit(const Loader& loader, const Gas& gas, Polygon* polygon)
{
    int index = polygon->getPhysicalIndex();
    const std::vector<std::string>& data = loader.getPhysicalData(index);
    polygon->f().f(gas.maxwell(data));
}


Integral IntegralConstructor(const Loader& loader)
{
    int intorder;
    try {
        intorder = loader.getData<int>("integral", "order");
    }
    catch (std::invalid_argument) {
        intorder = 1;
    }
    std::cout << "intorder = " << intorder << std::endl;

    bool is_free_molecular = !(loader.getData("integral", "value") == "True");

    int power = loader.getData<int>("integral", "power");

    SimpleSection* section;
    try {
        std::string str = loader.getData("integral", "section");
        if (str == "Simple") {
            throw std::invalid_argument("Simple section");
        }
        else if (str == "HS") {
            std::vector<double> ds;
            std::istringstream ss(loader.getData("integral", "ds"));
            double d;
            while (ss >> d)
                ds.push_back(d);
            section = new HSSection(ds);
            std::cout << "HSSection" << std::endl;

        }
        else if (str == "LJ") {
            std::vector<double> ds, es;
            std::istringstream ss(loader.getData("integral", "ds"));
            double x;
            while (ss >> x)
                ds.push_back(x);
            std::istringstream ss2(loader.getData("integral", "es"));
            while (ss2 >> x)
                es.push_back(x);
            std::string filename = loader.getData("integral", "file");
            section = new LJSection(ds, es, filename);
            std::cout << "LJSection" << std::endl;
        }
        else if (str == "Ma") {
            std::vector<double> as;
            std::istringstream ss(loader.getData("integral", "as"));
            double x;
            while (ss >> x)
                as.push_back(x);
            std::string filename = loader.getData("integral", "file");
            section = new MaSection(as, filename);
            std::cout << "MaSection" << std::endl;
        }
        else
            throw std::invalid_argument("Unknown section");
    }
    catch (std::invalid_argument) {
        section = new SimpleSection;
        std::cout << "SimpleSection" << std::endl;
    }

    return Integral(power, intorder, is_free_molecular, section);
}

*/

