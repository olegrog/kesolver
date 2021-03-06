#include <algorithm>
#include <stdexcept>
#include <sstream>

#include "Constructors.hpp"

#include "auxiliary.hpp"

#include "ximesh.hpp"
#include "ximesh_mixture.hpp"
#include "ximesh_mix.hpp"
#include "ximesh_rot.hpp"
#include "ximesh_rect.hpp"

#include "ci_simple.hpp"
#include "ci_mixture.hpp"
#include "ci_mix.hpp"
#include "ci_rect.hpp"

#include "base64/base64.hpp"
#include "logger/logger.hpp"

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
    std::cout << "masses.size() = " << masses.size() << std::endl;
    return masses;
}

double find_q(double d, int n, double err = 1e-10)
{
    // first approx
    if (d == n) return 1.;
    double q = 2 * d / n - 1;
    while (true) {
        double qq = std::pow(q, n-1);
        double f  = (q * qq - 1) / (q - 1);
        if (std::abs(f - d) < err)
            break;
        double df = ((n-1) * q * qq - n * qq + 1) / sqr(q-1);
        double dq = (d - f) / df;
        q += dq;
    }
    return q;
}

template <Symmetry symmetry>
std::pair<typename SymmetryTrait<symmetry>::VVd,
          typename SymmetryTrait<symmetry>::VVd>
makeVV(const typename SymmetryTrait<symmetry>::Vm v,
       const std::vector<double>& vs, 
       const std::vector<double>& vols,
       const std::vector<double>& vs2, 
       const std::vector<double>& vols2);

template <>
std::pair<SymmetryTrait<Cartesian>::VVd,
          SymmetryTrait<Cartesian>::VVd>
makeVV<Cartesian>(const SymmetryTrait<Cartesian>::Vm v,
                  const std::vector<double>& vs, 
                  const std::vector<double>& vols,
                  const std::vector<double>& vs2, 
                  const std::vector<double>& vols2)
{
    typedef SymmetryTrait<Cartesian>::VVd VVd;
    VVd vvs(vs2);
    for (int i = 0; i < 3; ++i)
        for (size_t j = 0; j < vvs[i].size(); ++j)
            vvs[i][j] += v[i];
    return std::make_pair(vvs, VVd(vols2));
}

template <>
std::pair<SymmetryTrait<Cylindrical>::VVd,
          SymmetryTrait<Cylindrical>::VVd>
makeVV<Cylindrical>(const SymmetryTrait<Cylindrical>::Vm v,
                    const std::vector<double>& vs, 
                    const std::vector<double>& vols,
                    const std::vector<double>& vs2, 
                    const std::vector<double>& vols2)
{
    typedef SymmetryTrait<Cylindrical>::VVd VVd;
    VVd vvs(vs2, vs);
    for (size_t j = 0; j < vvs[0].size(); ++j)
        vvs[0][j] += v;
    return std::make_pair(vvs, VVd(vols2, vols));
}

template <Symmetry symmetry>
std::pair<typename SymmetryTrait<symmetry>::VVd,
          typename SymmetryTrait<symmetry>::VVd>
makeVVByQ(const PropertyTree& tree, int rad, double cut,
          const typename SymmetryTrait<symmetry>::Vm v)
{
    double hcenter = tree["hcenter"].asDouble(); 
    double q = find_q(cut / hcenter, rad);
    std::cout << "q = " << q << std::endl;

    std::vector<double> vs, vols;
    double h = hcenter;
    double x = h / 2;
    for (int i = 0; i < rad; ++i) {
        vs.push_back(x);
        vols.push_back(h);
        x += h / 2;
        h *= q;
        x += h / 2;
    }

    std::vector<double> vs2, vols2;
    for (int i = rad-1; i >= 0; --i) {
        vs2.push_back(-vs[i]);
        vols2.push_back(vols[i]);
    }
    for (int i = 0; i < rad; ++i) {
        vs2.push_back(vs[i]);
        vols2.push_back(vols[i]);
    }

/*
    for (size_t i = 0; i < vs2.size(); ++i) 
        std::cout << vs2[i] << std::endl;
    for (size_t i = 0; i < vols2.size(); ++i) 
        std::cout << vols2[i] << std::endl;
*/

    return makeVV<symmetry>(v, vs, vols, vs2, vols2);
}

std::vector<double> readVectorFromBase64(const std::string& base64str)
{
    std::vector<unsigned char> bytes = base64::decode(base64str);
    const double* doubles = reinterpret_cast<const double*>(&bytes.front());
    const size_t size = bytes.size() / sizeof(double);
    std::vector<double> vector(size);
    for (size_t i = 0; i < size; ++i)
        vector[i] = doubles[i];
    return vector;
}

template <Symmetry symmetry>
typename SymmetryTrait<symmetry>::VVd
makeVFromRaw(int raw, const std::vector<double>& vs_raw);

template <>
SymmetryTrait<Cartesian>::VVd
makeVFromRaw<Cartesian>(int rad, const std::vector<double>& vs_raw)
{
    std::vector<double> vxs(2 * rad);
    std::vector<double> vys(2 * rad);
    std::vector<double> vzs(2 * rad);
    for (int i = 0; i < 2 * rad; ++i)
    {
        vxs[i] = vs_raw[i];
        vys[i] = vs_raw[i + 2 * rad];
        vzs[i] = vs_raw[i + 4 * rad];
    }
    return SymmetryTrait<Cartesian>::VVd(vxs, vys, vzs);
}

template <>
SymmetryTrait<Cylindrical>::VVd
makeVFromRaw<Cylindrical>(int rad, const std::vector<double>& vs_raw)
{
    std::vector<double> vxs(2 * rad);
    for (int i = 0; i < 2 * rad; ++i)
        vxs[i] = vs_raw[i];
    std::vector<double> vrs(rad);
    for (int i = 0; i < rad; ++i)
        vrs[i] = vs_raw[i + 2 * rad];
    return SymmetryTrait<Cylindrical>::VVd(vxs, vrs);
}

template <Symmetry symmetry>
std::pair<typename SymmetryTrait<symmetry>::VVd,
          typename SymmetryTrait<symmetry>::VVd>
makeVVread(const PropertyTree& tree, int rad)
{
    std::vector<double> vs_raw(readVectorFromBase64(tree["vs"].asString()));
    for (size_t i = 0; i < vs_raw.size(); ++i) 
        std::cout << vs_raw[i] << std::endl;
    typename SymmetryTrait<symmetry>::VVd vs(makeVFromRaw<symmetry>(rad, vs_raw));

    std::vector<double> vvs_raw(readVectorFromBase64(tree["vvs"].asString()));
    for (size_t i = 0; i < vvs_raw.size(); ++i) 
        std::cout << vvs_raw[i] << std::endl;
    typename SymmetryTrait<symmetry>::VVd vvs(makeVFromRaw<symmetry>(rad, vvs_raw));

    return std::make_pair(vs, vvs); 
}

template <Symmetry symmetry>
std::pair<typename SymmetryTrait<symmetry>::VVd,
          typename SymmetryTrait<symmetry>::VVd>
readVV(const PropertyTree& tree, int rad, double cut,
       const typename SymmetryTrait<symmetry>::Vm v)
{
    if (tree.isMember("hcenter"))
        return makeVVByQ<symmetry>(tree, rad, cut, v);
    else 
        return makeVVread<symmetry>(tree, rad);
}

template <Symmetry symmetry, Volume volume, TimeScheme scheme>
Gas* gasMixScheme(const PropertyTree& tree, 
                  const int rad, const double cut,
                  const typename SymmetryTrait<symmetry>::Vm v,
                  const std::vector<double>& masses)
{
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
                  const typename SymmetryTrait<symmetry>::Vm v,
                  const std::vector<double>& masses)
{
    std::string scheme;
    if (tree.isMember("time_scheme")) 
        scheme = tree["time_scheme"].asString();
    else
        scheme = "Continues";

    if      (scheme == "Euler") 
        return gasMixScheme<symmetry, volume, Euler>     (tree, rad, cut, v, masses);
    else if (scheme == "Continues")
        return gasMixScheme<symmetry, volume, Continues> (tree, rad, cut, v, masses);
    else
        throw std::invalid_argument("Unknown time scheme type");
}

template <Symmetry symmetry>
Gas* gasSymmetry(const PropertyTree& tree)
{
    LABEL

    const std::string type = tree["type"].asString();
    std::cout << "type = " << type << std::endl;

    int    rad = tree["rad"].asInt();
    double cut = tree["cut"].asDouble();
    std::cout << "rad, cut = " << rad << ' ' << cut << std::endl;

    typedef typename SymmetryTrait<symmetry>::Vm Vm;
    Vm v;
    if (tree.isMember("v")) 
        v = strTo<Vm>(tree["v"].asString());
    else 
        v = 0.;
    std::cout << "v = " << v << std::endl;

    if (type == "Simple") {
        LABEL

        std::string interp;
        if (tree.isMember("interpolation"))
            interp = tree["interpolation"].asString();
        else
            interp = "Power";

        if (interp == "Power") {
            typedef GasTemplate<symmetry, XiMesh, ColliderSimple<PowerInterp, symmetry> > GasType;
            return new GasType(XiMesh<symmetry>(rad, cut, v));
        } else if (interp == "PowerAndSymmetric") {
            typedef GasTemplate<symmetry, XiMesh, ColliderSimple<PowerAndSymmetricInterp, symmetry> > GasType;
            return new GasType(XiMesh<symmetry>(rad, cut, v));
        } else if (interp == "Linear") {
            typedef GasTemplate<symmetry, XiMesh, ColliderSimple<LinearInterp, symmetry> > GasType;
            return new GasType(XiMesh<symmetry>(rad, cut, v));
        } else if (interp == "Symmetric") {
            typedef GasTemplate<symmetry, XiMesh, ColliderSimple<SymmetricInterp, symmetry> > GasType;
            return new GasType(XiMesh<symmetry>(rad, cut, v));
        } else if (interp == "None") {
            typedef GasTemplate<symmetry, XiMesh, ColliderSimple<NoInterp, symmetry> > GasType;
            return new GasType(XiMesh<symmetry>(rad, cut, v));
        } else
            throw std::invalid_argument("Unknown interpolation type");

    }
    else if (type == "Mixture") {
        typedef GasTemplate<symmetry, XiMeshMixture, ColliderMixture<symmetry> > GasType;
        const std::vector<double> masses = readMasses(tree);

        const double a = cut / rad * std::sqrt(masses[0]);
        std::vector<int> rads;
        for (size_t i = 0; i < masses.size(); ++i) {
            std::cout << "mass[" << i << "] = " << masses[i] << std::endl;
            rads.push_back(static_cast<int>(rad * std::sqrt(masses[i] / masses[0])));
        }
        return new GasType(XiMeshMixture<symmetry>(a, rads, masses));
    }
    else if (type == "Mix") {
        const std::vector<double> masses = readMasses(tree);
        const std::string volume = tree["volume"].asString();
        if      (volume == "Wide") 
            return gasMixVolume<symmetry, Wide>      (tree, rad, cut, v, masses);
        else if (volume == "Symmetric")
            return gasMixVolume<symmetry, Symmetric> (tree, rad, cut, v, masses);
        else if (volume == "Tight")
            return gasMixVolume<symmetry, Tight>     (tree, rad, cut, v, masses);
        else if (volume == "Tight2")
            return gasMixVolume<symmetry, Tight2>    (tree, rad, cut, v, masses);
        else if (volume == "Grad")
            return gasMixVolume<symmetry, Grad>      (tree, rad, cut, v, masses);
        else
            throw std::invalid_argument("Unknown mix volume type");
    }
    else if (type == "Rect") {
        typedef typename SymmetryTrait<symmetry>::VVd VVd;
        std::pair<VVd, VVd> vpair = readVV<symmetry>(tree, rad, cut, v);
        VVd vv(vpair.first), vvol(vpair.second);

        XiMeshRect<symmetry> ximesh(cut, v, vv, vvol);

        const std::string volume = 
            tree.isMember("volume") ? tree["volume"].asString() : "Tight";

        if (volume == "Symmetric")
            return new GasTemplate <
                                     symmetry, XiMeshRect,
                                     ColliderRect<symmetry, Symmetric>
                                   > (ximesh);
        else if (volume == "Tight")
            return new GasTemplate <
                                     symmetry, XiMeshRect,
                                     ColliderRect<symmetry, Tight>
                                   > (ximesh);
//        else if (volume == "Grad")
//            return new GasTemplate <
//                                     symmetry, XiMeshRect, 
//                                     ColliderRect<symmetry, Grad> 
//                                   > (ximesh);
        else
            throw std::invalid_argument("Unknown rect volume type");
    }
    else
        throw std::invalid_argument("Unknown gas type");
}

void GasConstructor(const PropertyTree& tree, Gas** gas_pp)
{
    LABEL

    std::string symm;
    if (tree.isMember("symmetry"))
        symm = tree["symmetry"].asString();
    else 
        symm = "Cartesian";
    std::cout << "symm = " << symm << std::endl;

    Gas* gas_p;
    if (symm == "Cylindrical")
        gas_p = gasSymmetry<Cylindrical>(tree);
    else if (symm == "Cartesian")
        gas_p = gasSymmetry<Cartesian>(tree);
    else
        throw std::invalid_argument("Unknown symmetry type");

    *gas_pp = gas_p;
}

void GivePolygonMemoryAndInit(const PropertyTree& tree, const Gas& gas, Polygon* polygon)
{
    const std::string name   = polygon->getPhysicalName();
    const PropertyTree& data = tree["initial_conditions"][name];
    const std::string type   = data["type"].asString();
    const std::string index  = std::to_string(polygon->getIndex());
    if (type == "maxwell") {
        polygon->f().f(gas.maxwell(data));
    }
    else if (type == "maxwell-nonuniform") {
        polygon->f().f(gas.maxwell(data["values"][index]));
    }
    else if (type == "grad13-nonuniform") {
        polygon->f().f(gas.grad13(data["values"][index]));
    }
    else if (type == "raw") {
        std::vector<unsigned char> bytes = base64::decode(data["raw"].asString());
        const double* doubles = reinterpret_cast<const double*>(&bytes.front());
        const size_t size = bytes.size() / sizeof(double);
        LOG(INFO) << "bytes.size() = " << bytes.size();
        LOG(INFO) << "f().size() = " << size;
        std::vector<double> f(size);
        for (size_t i = 0; i < size; ++i)
            f[i] = doubles[i];
        polygon->f().f(f);
    }
    else 
        throw std::invalid_argument("Unknown initial condition type");
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
            while (ss >> d) {
                ds.push_back(d);
                std::cout << "d = " << d << std::endl;
            }
            section = new HSSection(ds);
            std::cout << "HSSection" << std::endl;
        }
        else if (str == "VHS") {
            std::vector<double> ds;
            std::istringstream ss(tree["ds"].asString());
            double d;
            while (ss >> d) {
                ds.push_back(d);
                std::cout << "d = " << d << std::endl;
            }
            double omega = tree["omega"].asDouble();
            section = new VHSSection(ds, omega);
            std::cout << "VHSSection" << std::endl;

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
        else if (str == "Anikin") {
            double d = strTo<double>(tree["d"].asString());
            double e = strTo<double>(tree["e"].asString());

            std::string file_teta  = tree["file_teta"].asString();
            std::string file_vel   = tree["file_vel"].asString();
            std::string file_sigma = tree["file_sigma"].asString();
            section = new AnikinSection(d, e, file_teta, file_vel, file_sigma);
            std::cout << "AnikinSection" << std::endl;
        }
        else if (str == "AbInitio") {
            double d    = tree["d"].asDouble();
            double temp = tree["temp"].asDouble();

            std::vector< std::string > filenames;
            std::istringstream ss(tree["file"].asString());
            std::string x;
            while (ss >> x)
                filenames.push_back(x);

            section = new AbInitioSection(d, temp, filenames);
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

