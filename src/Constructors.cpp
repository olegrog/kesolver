#include <algorithm>
#include <stdexcept>
#include "Loader.hpp"
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

Polygon* createPolygon(const std::vector<int>& data) {
    if      (data[0] == 4) return new Tetrahedron();
    else if (data[0] == 5) return new Hexahedron();
    else if (data[0] == 6) return new Prism();
    else return 0;
}

void constructPolygon(const std::vector<int>& data, const std::vector<V3d>& nodes, 
        const std::vector<Polygon*>& polygons, Polygon* polygon) {

    std::vector<V3d> vertexes(polygon->getNumberOfVertexes());

    std::vector<int>::const_iterator data_p = ++data.begin();
    for (size_t i = 0; i < vertexes.size(); ++i) 
        vertexes[i] = nodes[*data_p++];

    std::vector<int> neigbors(*data_p++);
    for (size_t i = 0; i < neigbors.size(); ++i) 
        neigbors[i] = *data_p++;
    
    int rank = *data_p++;
    int physical_index = *data_p++;

    polygon->setVertexCoordinates(vertexes);
    polygon->calculateLength();
    polygon->calculateVolume();
    polygon->calculateCenter();
//  std::cout << polygon->getVolume() << ' ' << polygon->getCenter() << std::endl;
    polygon->setNeigbors(neigbors);
    polygon->setRank(rank);
    polygon->setPhysicalIndex(physical_index);

}

void PolygonsConstructor(const Loader& loader, std::vector<Polygon*>& polygons) {

    for (std::vector< std::vector<int> >::const_iterator p = loader.getCells().begin();
            p != loader.getCells().end(); ++p) 
        polygons.push_back(createPolygon(*p));
    
    std::vector< std::vector<int> >::const_iterator data_p = loader.getCells().begin();
    for (std::vector<Polygon*>::iterator p = polygons.begin(); p != polygons.end(); ++p) 
        constructPolygon(*data_p++, loader.getNodes(), polygons, *p);
}

PhysicalFacet* createFacet(const std::vector<std::string>& data, const Gas& gas) {
    if (data[0] == "w")  {
        std::vector<std::string> newdata;
        std::copy(++data.begin(), data.end(), std::back_inserter(newdata));
        return new WallMaxwellFacet(gas.maxwell(newdata));
    }
    else if (data[0] == "l")  {
        std::vector<std::string> newdata;
        std::copy(++data.begin(), data.end(), std::back_inserter(newdata));
        return new MaxwellFacet(gas.maxwell(newdata));
    }
    else if (data[0] == "m") {
        return new MirrorFacet((Axis)strTo<int>(data[1]));
    }
    else {
        return new GateFacet();
    }
}

bool lessFacet(PhysicalFacet* f1, PhysicalFacet* f2) {
    return f1->order() < f2->order();
}

void FacetsConstructor(const Loader& loader, 
                       const std::vector<Polygon*>& spacemesh, 
                       std::vector<PhysicalFacet*>& facets, 
                       const Gas& gas,
                       double time_step, int rank)
{
    for (std::vector< std::vector<int> >::const_iterator pp = loader.getFacets().begin();
            pp != loader.getFacets().end(); ++pp) {

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
        facet->findNormalAndSquare();
//      std::cout << facet->getSquare() << ' ' << facet->getCenter() << std::endl;
        facet->findMultInOut(time_step, spacemesh);

        facets.push_back(facet);
    }
    std::sort(facets.begin(), facets.end(), lessFacet);
}

void MypolysConstructor(int rank, const std::vector<Polygon*>& polygons, std::vector<int>& mypolys) {
    for (size_t i = 0; i < polygons.size(); ++i) 
        if (polygons[i]->getRank() == rank) {
            mypolys.push_back(i);
        }
}


void GivePolygonMemoryAndInit(const Loader& loader, const Gas& gas, Polygon* polygon)
{
    int index = polygon->getPhysicalIndex();
    const std::vector<std::string>& data = loader.getPhysicalData(index);
    polygon->f().f(gas.maxwell(data));
}

double findTimeStep(const std::vector<Polygon*>& spacemesh, const Gas& gas, double curnt) {
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

const std::vector<double> readMasses(const Loader& loader) {
    std::vector<double> masses;
    try {
        const double m1   = 1;
        const double m2   = loader.getData<double>("gas", "mass_ratio");
        const double m0   = 0.5 * (m1 + m2);
        masses.push_back(m1 / m0);
        masses.push_back(m2 / m0);
    }
    catch (std::invalid_argument) {
        std::istringstream ss(loader.getData("gas", "masses"));
        double mass;
        while (ss >> mass)
            masses.push_back(mass);
    }
    return masses;
}

template <Symmetry symmetry, Volume volume, TimeScheme scheme>
Gas* gasMixScheme(const Loader& loader, 
                  const int rad, const double cut,
                  const std::vector<double>& masses)
{
    typedef typename SymmetryTrait<symmetry>::Vm Vm;
    Vm v;
    try {
        v = loader.getData<Vm>("gas", "v");
    }
    catch (std::invalid_argument) {
        v = 0.;
    }
    std::vector<double> adds;
    try {
        std::istringstream ss(loader.getData("gas", "adds"));
        double add;
        while (ss >> add)
            adds.push_back(add);
    }
    catch (std::invalid_argument) {
        for (size_t i = 0; i < masses.size(); ++i)
            adds.push_back(0.0);
    }
    typedef GasTemplate<symmetry, XiMeshMix, ColliderMix<symmetry, volume, scheme> > GasType;
    return new GasType(XiMeshMix<symmetry>(rad, cut, masses, v, adds));
}

template <Symmetry symmetry, Volume volume>
Gas* gasMixVolume(const Loader& loader, 
                  const int rad, const double cut,
                  const std::vector<double>& masses)
{
    std::string scheme;
    try {
        scheme = loader.getData("gas", "time_scheme");
    }
    catch (std::invalid_argument) {
        scheme = "Continues";
    }

    if      (scheme == "Euler") 
        return gasMixScheme<symmetry, volume, Euler>     (loader, rad, cut, masses);
    else if (scheme == "Continues")
        return gasMixScheme<symmetry, volume, Continues> (loader, rad, cut, masses);
    else { return NULL; /* TODO raise */ }
}

template <Symmetry symmetry>
Gas* gasSymmetry(const Loader& loader, const std::string& type, 
                 const int rad, const double cut)
{
    LABEL
    if      (type == "Simple") {
        typedef typename SymmetryTrait<symmetry>::Vm Vm;
        Vm v;
        try {
            v = loader.getData<Vm>("gas", "v");
        }
        catch (std::invalid_argument) {
            v = 0.;
        }
        typedef GasTemplate<symmetry, XiMesh, ColliderSimple<symmetry> > GasType;
        return new GasType(XiMesh<symmetry>(rad, cut, v));
    }
    else if (type == "Mixture") {
        typedef GasTemplate<symmetry, XiMeshMixture, ColliderMixture<symmetry> > GasType;
        const std::vector<double> masses = readMasses(loader);

        const double a = cut / rad * std::sqrt(masses[0]);
        std::vector<int> rads;
        for (size_t i = 0; i < masses.size(); ++i)
            rads.push_back(static_cast<int>(rad * std::sqrt(masses[i] / masses[0])));
        return new GasType(XiMeshMixture<symmetry>(a, rads, masses));
    }
    else if (type == "Mix") {
        const std::vector<double> masses = readMasses(loader);
        const std::string volume = loader.getData("gas", "volume");
        if      (volume == "Wide") 
            return gasMixVolume<symmetry, Wide>      (loader, rad, cut, masses);
        else if (volume == "Symmetric")
            return gasMixVolume<symmetry, Symmetric> (loader, rad, cut, masses);
        else if (volume == "Tight")
            return gasMixVolume<symmetry, Tight>     (loader, rad, cut, masses);
        else if (volume == "Tight2")
            return gasMixVolume<symmetry, Tight2>    (loader, rad, cut, masses);
        else if (volume == "Grad")
            return gasMixVolume<symmetry, Grad>      (loader, rad, cut, masses);
        else { return NULL; /* TODO raise */ }
    }
    else { return NULL; /* TODO raise */ }
}

void GasConstructor(const Loader& loader, Gas** gas_pp)
{
    const std::string type = loader.getData("gas", "type");
    std::cout << "type = " << type << std::endl;

    int    rad = loader.getData<int>("gas", "rad");
    double cut = loader.getData<double>("gas", "cut");
    std::cout << "rad, cut = " << rad << ' ' << cut << std::endl;

    std::string symm;
    try {
        symm = loader.getData("gas", "symmetry");
    }
    catch (std::invalid_argument) {
        symm = "Cartesian";
    }
    std::cout << "symm = " << symm << std::endl;

    Gas* gas_p;
    if      (symm == "Cylindrical")
        gas_p = gasSymmetry<Cylindrical>(loader, type, rad, cut);
    else if (symm == "Cartesian")
        gas_p = gasSymmetry<Cartesian>(loader, type, rad, cut);
    else { gas_p = NULL; /* TODO raise */ }

    *gas_pp = gas_p;
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

