#include <algorithm>
#include <stdexcept>
#include <sstream>

#include "Constructors.hpp"

#include "auxiliary.hpp"

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



void GivePolygonMemoryAndInit(const PropertyTree& tree, const Gas& gas, Polygon* polygon)
{
    const std::string name = polygon->getPhysicalName();
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

