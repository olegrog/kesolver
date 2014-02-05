#include <iostream>
#include <fstream>

#include "section.hpp"

LJSection::LJSection(const std::vector<double>& ds,
                     const std::vector<double>& es,
                     const std::string& filename) : ds(ds), es(es)
{
    std::ifstream fd(filename.c_str());
    size_t m, n;
    fd >> m >> n;
    gs.reserve(m);
    ths.reserve(n);
    sigma.reserve(m*n);
    for (size_t i = 0; i < m; ++i) {
        double x;
        fd >> x;
        gs.push_back(x);
    }
    for (size_t i = 0; i < n; ++i) {
        double x;
        fd >> x;
        ths.push_back(x);
    }
    for (size_t i = 0; i < m; ++i) 
        for (size_t j = 0; j < n; ++j) {
            double x;
            fd >> x;
            sigma.push_back(x);
        }
    fd.close();
}

int bin_search(const double x, const std::vector<double>& xs)
{
    if (x < xs[0]) return 0;
    if (xs[xs.size()-1] < x) return xs.size() - 1;

    size_t l = 0;
    size_t r = xs.size() - 1;

    size_t i = 0;
    while (true) {
        i = (l + r) / 2;
        if (xs[i] > x)
            r = i;
        else if (xs[i+1] < x)
            l = i + 1;
        else
            break;
    }

    if (x - xs[i] < xs[i+1] - x) return i;
    else return i + 1;
}

double LJSection::ker(const double g, const double cth) const
{
    int i = bin_search(g, gs);
    double th = std::acos(cth);
    int j = bin_search(th, ths);
/*
    std::cout << "g, i = " << g << ' ' << i << std::endl;
    std::cout << "th, j = " << th << ' ' << j << std::endl;
    std::cout << "sigma = " << sigma[i * ths.size() + j] << std::endl;
*/
    return sigma[i * ths.size() + j]; // / M_PI;
}

MaSection::MaSection(const std::vector<double>& as,
                     const std::string& filename) : as(as)
{
    std::ifstream fd(filename.c_str());
    size_t n;
    fd >> n;
    ths.reserve(n);
    sigma.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        double x;
        fd >> x;
        ths.push_back(x);
    }
    for (size_t j = 0; j < n; ++j) {
        double x;
        fd >> x;
        sigma.push_back(x);
    }
    fd.close();
}


double MaSection::ker(const double g, const double cth) const
{
    if (std::abs(g) < 1e-12)
        return 0;
    double th = std::acos(cth);
    int j = bin_search(th, ths);
//    std::cout << "th, g, sigma = " << th << ' ' << g << ' ' << sigma[j] << ' ' << sigma[j] / M_PI / g << std::endl;
    return sigma[j] / M_PI / g;
}

AnikinSection::AnikinSection(double d, double e,
                             const std::string& file_teta,
                             const std::string& file_vel,
                             const std::string& file_sigma)
    : d(d), e(e)
{
    std::ifstream fd_teta(file_teta.c_str());
    double x;
    while (fd_teta >> x)
        ths.push_back(x);

    std::ifstream fd_vel(file_vel.c_str());
    while (fd_vel >> x)
        gs.push_back(x);

    std::ifstream fd_sigma(file_sigma.c_str());
    while (fd_sigma >> x)
        sigma.push_back(x);
}

double AnikinSection::section(const double g, const double cth) const 
{
    return section(g, cth, 0, 0);
}

double AnikinSection::section(const double g, const double cth,
                              const int i1, const int i2) const
{
    double g1 = g / std::sqrt(e);

    int i = bin_search(g1, gs);
    double th = std::acos(cth);
    int j = bin_search(th, ths);    

    int l;
    if (i1 == i2) l = i1;
    else          l = 2 + i1 + i2; // TODO 2

    return sigma[l * ths.size() * gs.size() + i * ths.size() + j]
        / std::sin(th) * sqr(d);
}

AbInitioSection::AbInitioSection(double diam,
                     double temp,
                     const std::vector<std::string>& sigmafilenames)
: diam_(diam), temp_(temp)
{
    size_t number_of_crosses = 0;
    int n = 0;
    while (number_of_crosses < sigmafilenames.size()) {
        n++;
        number_of_crosses += n;
    }
    if (number_of_crosses != sigmafilenames.size())
       ; // TODO raise exception

    number_of_components_ = n;

    sigmas.resize(number_of_crosses);
    for (size_t i = 0; i < number_of_crosses; ++i) 
    {
        const std::string& filename = sigmafilenames[i];
        std::vector<double>& gs     = sigmas[i].gs;
        std::vector<double>& ths    = sigmas[i].ths;
        std::vector<double>& sigma  = sigmas[i].sigma;

        std::ifstream fd(filename.c_str());
        size_t m, n;
        fd >> m >> n;

        gs.reserve(m);
        ths.reserve(n);
        sigma.reserve(m*n);
        for (size_t i = 0; i < m; ++i) {
            double x;
            fd >> x;
            gs.push_back(x);
        }
        for (size_t i = 0; i < n; ++i) {
            double x;
            fd >> x;
            ths.push_back(x);
        }
        for (size_t i = 0; i < m; ++i) 
            for (size_t j = 0; j < n; ++j) {
                double x;
                fd >> x;
                sigma.push_back(x);
            }
        fd.close();
    }

}

double AbInitioSection::ker(const double g, const double cth, const int i1, const int i2) const
{
    const int j1 = std::min(i1, i2);
    const int j2 = std::max(i1, i2);

    int k = j2 - j1 + (number_of_components_ * j1 - (j1 - 1) * j1 / 2);

    double g1 = g * sqrt(temp_);
    int i = bin_search(g1, sigmas[k].gs);
    double th = std::acos(cth);
    int j = bin_search(th, sigmas[k].ths);

    double sigma = sigmas[k].sigma[i * sigmas[k].ths.size() + j];

    return  sigma * sqr(diam_);
}

