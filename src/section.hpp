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
        virtual double ker(const double g, const double th) const { return 0.25; }
        std::vector<double> ds;
};

class VHSSection : public HSSection {
    public:
        VHSSection(const std::vector<double>& ds, double omega) : HSSection(ds), omega(omega) {}
    private:
        virtual double ker(const double g, const double th) const { return 0.25*std::pow(g, .5-omega); }
        double omega;
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

class AnikinSection : public SimpleSection {
    public:
        AnikinSection(double d, double e, 
                      const std::string& file_teta,
                      const std::string& file_vel,
                      const std::string& file_sigma);

        virtual double section(const double g, const double th) const;
        virtual double section(const double g, const double th,
                               const int i1, const int i2) const;

    private:
        double d, e;
        std::vector<double> gs, ths;
        std::vector<double> sigma;
};

class AbInitioSection : public SimpleSection {
    public:
        AbInitioSection(const double diam,
                  const double temp,
                  const std::vector<std::string>& sigmafilenames);
        virtual double section(const double g, const double th) const {
            return ker(g, th, 0, 0);
        }
        virtual double section(const double g, const double th,
                               const int i1, const int i2) const {
            return ker(g, th, i1, i2);
        }
    private:
        double ker(const double g, const double th, int i1, int i2) const;

        double diam_, temp_;
        int number_of_components_;

        struct SigmaData {
            std::vector<double> gs, ths;
            std::vector<double> sigma;
        };

        std::vector<SigmaData> sigmas;
};

#endif
