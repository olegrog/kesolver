#ifndef _MACRO_HPP_
#define _MACRO_HPP_

#include <boost/tuple/tuple.hpp>
#include <string>

#include "ximesh.hpp"
#include "ximesh_mixture.hpp"
#include "ximesh_mix.hpp"
#include "ximesh_rot.hpp"
#include "ximesh_rect.hpp"

template <Symmetry symmetry>
struct MacroSimple {
    typedef typename SymmetryTrait<symmetry>::Vm Vm;
    typedef typename SymmetryTrait<symmetry>::Vd Vd;
    double n;
    Vm v;
    double temp;
    Vd t;
    Vm q;
    double h;

    MacroSimple(const double n, const Vm v,
                const double temp, const Vd t, 
                const Vm q,
                const double h) :
            n(n), v(v),
            temp(temp), t(t),
            q(q), h(h) {}

    template <typename F, template <Symmetry > class XiMeshType>
    MacroSimple(const F& f, const XiMeshType<symmetry>& mesh)
    {
        n = 0;
        temp = 0;
        v = 0;
        q = 0;
        h = 0;
        for (size_t i = 0; i < mesh.size(); ++i) {
            n    += f[i];
            v    += mesh.p(i) * f[i];
            temp += mesh.e(i) * f[i];
            t    += mesh[i] * mesh[i] * f[i];
            q    += 0.5 * mesh.p(i) * mesh.e(i) * f[i];
            h    += f[i] * std::log(f[i]);
        }
        v /= n;
        temp = (temp - n * sqr(v)) / 3 / n;
        Vd u = vm2vd(v);
        t    = (t    - n * u * u) * SymmetryTrait<symmetry>::t_mult() / n;
        n *= mesh.vol();
        q *= mesh.vol();
        h *= mesh.vol();
    }

};

template <Symmetry symmetry>
struct MacroMixture {
    typedef typename MacroSimple<symmetry>::Vm Vm;
    typedef typename MacroSimple<symmetry>::Vd Vd;
    double n;
    Vm v;
    double temp;
    Vd t;
    Vm q;
    double h;
    std::vector< MacroSimple<symmetry> > macros;

    template <typename F, template <Symmetry > class XiMeshType>
    MacroMixture(const F& f, const XiMeshType<symmetry>& mesh)
    {
        boost::tie(n, v, temp, t, q, h) = macro(f, mesh, 0, mesh.size());

        macros.reserve(mesh.xiMeshes().size());
        for (size_t j = 0; j < mesh.xiMeshes().size(); ++j) {
            double n, temp;
            Vm v;
            Vd t;
            Vm q;
            double h;
            boost::tie(n, v, temp, t, q, h) = macro(f, mesh,
                                         mesh.offset(j), mesh.offset(j+1));
            macros.push_back(MacroSimple<symmetry>(n, v, temp, t, q, h));
        }
    }

    private:
        
        template <typename F, template <Symmetry > class XiMeshType>
        boost::tuple<double, Vm, double, Vd, Vm, double> macro(const F&  f,
                                             const XiMeshType<symmetry>& mesh,
                                             const int i1, const int i2)
        {
            double n = 0, m = 0, temp = 0;
            Vm v = 0;
            Vd t = 0;
            Vm q = 0;
            double h = 0;
            for (int i = i1; i < i2; ++i) {
                n    += f[i];
                m    += f[i] * mesh.m(i);
                v    += f[i] * mesh.p(i);
                temp += f[i] * sqr(mesh[i]);
                t    += f[i] * mesh.m(i) * mesh[i] * mesh[i];
                q    += f[i] * 0.5 * mesh.p(i) * mesh.e(i);
                h    += f[i] * log(f[i]);
            }
            v /= m;
            q *= n / m;
            temp = (temp - m * sqr(v)) / 3 / n;
            Vd u = vm2vd(v);
            t    = (t    - m * u * u) * SymmetryTrait<symmetry>::t_mult() / n;
            n *= mesh.vol();
            q *= mesh.vol();
            h *= mesh.vol();
            return boost::make_tuple(n, v, temp, t, q, h);
        }

};

template <Symmetry symmetry>
struct MacroRot {
    typedef typename SymmetryTrait<symmetry>::Vm Vm;
    typedef typename SymmetryTrait<symmetry>::Vd Vd;
    double n;
    Vm v;
    double temp, t_rot;
    Vd t;
    Vm q;
    double h;

    template <typename F, template <Symmetry > class XiMeshType>
    MacroRot(const F& f, const XiMeshType<symmetry>& mesh)
    {
        n = 0;
        temp = 0; t_rot = 0; t = 0;
        v = 0;
        q = 0;
        h = 0;
        for (size_t i = 0; i < mesh.size(); ++i) {
            n     += f[i];
            v     += mesh.p(i) * f[i];
            temp  += mesh.e(i) * f[i];
            t     += mesh[i] * mesh[i] * f[i];
            t_rot += mesh.erot(i) * f[i];
            q     += mesh.p(i) * mesh.e(i) * f[i];
            h     += f[i] * std::log(f[i]);
        }
        v     /= n;
        temp   = (temp - n * sqr(v)) / 5 / n;
        t_rot /= 2 * n;
        Vd u   = vm2vd(v);
        t      = (t    - n * u * u) * SymmetryTrait<symmetry>::t_mult() / n;
        n     *= mesh.vol();
        q     *= mesh.vol();
        h     *= mesh.vol();
    }
        
};

template <Symmetry symmetry> inline
const std::string toStr(const MacroSimple<symmetry>& s) {
    std::ostringstream ss;
    ss << "[ " << s.n << ' ' << s.v << ' ' << s.temp << ' ' << s.t << ' ' << s.q << ' ' << s.h << " ]";
    return ss.str();
}

template <Symmetry symmetry> inline
const std::string toStr(const MacroMixture<symmetry>& s) {
    std::ostringstream ss;
    ss <<  "[ " << 
            s.n << ' ' << s.v << ' ' << 
            s.temp << ' ' << s.t << ' ' << s.q << ' ' << s.h << ' ';
        for (size_t i = 0; i < s.macros.size(); ++i)
            ss << ' ' << s.macros[i];
    ss << " ]";
    return ss.str();
}
        
template <Symmetry symmetry> inline
const std::string toStr(const MacroRot<symmetry>& s) {
    std::ostringstream ss;
    ss << "[ " << 
           s.n << ' ' << s.v << ' ' << 
           s.temp << ' ' << s.t << ' ' << s.t_rot << ' ' << 
           s.q << ' ' << s.h << 
         " ]";
    return ss.str();
}


template <Symmetry symmetry> inline
std::ostream& operator<<(std::ostream& to, const MacroSimple<symmetry>& s)
{
    to << toStr(s);
    return to;
}

template <Symmetry symmetry> inline
std::ostream& operator<<(std::ostream& to, const MacroMixture<symmetry>& s)
{
    to << toStr(s);
    return to;
}

template <Symmetry symmetry> inline
std::ostream& operator<<(std::ostream& to, const MacroRot<symmetry>& s)
{
    to << toStr(s);
    return to;
}

template <typename F, Symmetry symmetry>
MacroSimple<symmetry> macro(const F& f, const XiMesh<symmetry>& mesh)
{
    return MacroSimple<symmetry>(f, mesh);
}

template <typename F, Symmetry symmetry>
MacroMixture<symmetry> macro(const F& f, const XiMeshMixture<symmetry>& ximesh)
{
    return MacroMixture<symmetry>(f, ximesh);
}

template <typename F, Symmetry symmetry>
MacroMixture<symmetry> macro(const F& f, const XiMeshMix<symmetry>& ximesh)
{
    return MacroMixture<symmetry>(f, ximesh);
}

template <typename F, Symmetry symmetry>
MacroRot<symmetry> macro(const F& f, const XiMeshRot<symmetry>& mesh)
{
    return MacroRot<symmetry>(f, mesh);
}

template <typename F, Symmetry symmetry>
MacroSimple<symmetry> macro(const F& f, const XiMeshRect<symmetry>& mesh)
{
    return MacroSimple<symmetry>(f, mesh);
}

template <typename DF1, typename DF2, typename V, typename XiMeshType>
void equateStreamsSimple(DF1& g, const DF2& f, const V n, const XiMeshType& mesh)
{
    double sf = 0., sg = 0.;
    for (size_t i = 0; i < mesh.size(); i++) {
        double xin = dot(mesh[i], n);
        if (xin < 0.0) {
            sf += f[i] * xin;
            sg += g[i] * xin;
        }
    }
    for (size_t i = 0; i < mesh.size(); i++)
        g[i] *= sf / sg;
}

template <typename DF, typename V, typename XiMeshType>
void equateStreamsMixture(DF& g, const DF& f, const V n, const XiMeshType& mesh)
{
    static std::vector<double> sf, sg;
    size_t size = mesh.xiMeshes().size();
    sf.resize(size);
    sg.resize(size);
    for (size_t i = 0; i < size; ++i) {
        sf[i] = 0.;
        sg[i] = 0.;
    }

    for (size_t i = 0; i < mesh.size(); i++) {
        double xin = dot(mesh[i], n);
        if (xin < 0.0) {
            sf[mesh.i2ci(i)] += f[i] * xin;
            sg[mesh.i2ci(i)] += g[i] * xin;
        }
    }
/*
    for (size_t i = 0; i < size; ++i) {
        std::cout << "s[" << i << "] = " << sf[i] << ' ' << sg[i] << std::endl;
    }
*/
    for (size_t i = 0; i < mesh.size(); i++)
        g[i] *= sf[mesh.i2ci(i)] / sg[mesh.i2ci(i)];

}

template <typename DF, typename V, Symmetry symmetry>
void equateStreams(DF& g, const DF& f, const V n, const XiMesh<symmetry>& mesh)
{
    return equateStreamsSimple(g, f, n, mesh);
}

template <typename DF, typename V, Symmetry symmetry>
void equateStreams(DF& g, const DF& f, const V n, const XiMeshMixture<symmetry>& mesh)
{
    return equateStreamsMixture(g, f, n, mesh);
}

template <typename DF, typename V, Symmetry symmetry>
void equateStreams(DF& g, const DF& f, const V n, const XiMeshMix<symmetry>& mesh)
{
    return equateStreamsMixture(g, f, n, mesh);
}

template <typename DF, typename V, Symmetry symmetry>
void equateStreams(DF& g, const DF& f, const V n, const XiMeshRot<symmetry>& mesh)
{
    return equateStreamsSimple(g, f, n, mesh);
}

template <typename DF, typename V, Symmetry symmetry>
void equateStreams(DF& g, const DF& f, const V n, const XiMeshRect<symmetry>& mesh)
{
    return equateStreamsSimple(g, f, n, mesh);
}

template <typename DF, typename Stream, typename XiMeshType>
Stream& equateNSimple(DF& g, Stream& stream, const XiMeshType& mesh)
{
    double n;
    stream >> n;
    double z = 0.;
    for (size_t i = 0; i < mesh.size(); ++i)
        z += g[i];
    z *= mesh.vol() / n;
    for (size_t i = 0; i < mesh.size(); ++i) 
        g[i] /= z;
    return stream;
}

template <typename DF, typename Stream, typename XiMeshType>
Stream& equateNMixture(DF& g, Stream& stream, const XiMeshType& mesh)
{
    static std::vector<double> zs;
    static std::vector<double> ns;
    size_t size = mesh.xiMeshes().size();
    zs.resize(size);
    ns.resize(size);
    for (size_t i = 0; i < size; ++i) {
        zs[i] = 0.;
        ns[i] = 1.;
        stream >> ns[i];
    }
    for (size_t i = 0; i < mesh.size(); i++)
        zs[mesh.i2ci(i)] += g[i];
    for (size_t i = 0; i < size; ++i)
        zs[i] *= mesh.vol() / ns[i];
    for (size_t i = 0; i < mesh.size(); i++)
        g[i] /= zs[mesh.i2ci(i)];
    return stream;
}

template <typename DF, typename Stream, Symmetry symmetry>
Stream& equateN(DF& g, Stream& stream, const XiMesh<symmetry>& mesh)
{
    return equateNSimple(g, stream, mesh);
}

template <typename DF, typename Stream, Symmetry symmetry>
Stream& equateN(DF& g, Stream& stream, const XiMeshMixture<symmetry>& mesh)
{
    return equateNMixture(g, stream, mesh);
}

template <typename DF, typename Stream, Symmetry symmetry>
Stream& equateN(DF& g, Stream& stream, const XiMeshMix<symmetry>& mesh)
{
    return equateNMixture(g, stream, mesh);
}

template <typename DF, typename Stream, Symmetry symmetry>
Stream& equateN(DF& g, Stream& stream, const XiMeshRot<symmetry>& mesh)
{
    return equateNSimple(g, stream, mesh);
}

template <typename DF, typename Stream, Symmetry symmetry>
Stream& equateN(DF& g, Stream& stream, const XiMeshRect<symmetry>& mesh)
{
    return equateNSimple(g, stream, mesh);
}

#endif
