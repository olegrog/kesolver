#ifndef _SECTION_HPP_
#define _SECTION_HPP_

#include <vector>
#include <string>

#include "auxiliary.hpp"

class SimpleSection {
    public:
        virtual double section(const double g, const double th) const { return 0.25; }
        virtual double section(const double g, const double th,
                               const int i1, const int i2) const {
            return section(g, th);
        }
};

class HSSection : public SimpleSection {
    public:
        HSSection(const std::vector<double>& ds) : ds(ds) {}
        virtual double section(const double g, const double th) const {
            return ker(g, th) * sqr(ds[0]);
        }
        virtual double section(const double g, const double th,
                               const int i1, const int i2) const {
            return ker(g, th) * sqr(ds[i1] + ds[i2]) / 4;
        }
    private:
        double ker(const double g, const double th) const { return 0.25; }
        std::vector<double> ds;
};

class LJSection : public SimpleSection {
    public:
        LJSection(const std::vector<double>& ds,
                  const std::vector<double>& es,
                  const std::string& filename);
        virtual double section(const double g, const double th) const {
            double g1 = g / std::sqrt( es[0] );
            return ker(g1, th) * sqr(ds[0]);
        }
        virtual double section(const double g, const double th,
                               const int i1, const int i2) const {
            double g1 = g / std::sqrt( std::sqrt( es[i1] * es[i2] ) );
            return ker(g1, th) * sqr(ds[i1] + ds[i2]) / 4;
        }
    private:
        double ker(const double g, const double th) const;

        std::vector<double> ds;
        std::vector<double> es;

        std::vector<double> gs, ths;
        std::vector<double> sigma;
};

class MaSection : public SimpleSection {
    public:
        MaSection(const std::vector<double>& as,
                  const std::string& filename);
        virtual double section(const double g, const double th) const {
            double g1 = g / std::sqrt( as[0] );
            return ker(g1, th);
        }
        virtual double section(const double g, const double th,
                               const int i1, const int i2) const {
            double g1 = g / std::sqrt( std::sqrt( as[i1] * as[i2] ) );
            return ker(g1, th);
        }
    private:
        double ker(const double g, const double th) const;

        std::vector<double> as;

        std::vector<double> ths;
        std::vector<double> sigma;
};

#endif
