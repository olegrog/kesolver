#include <exception>
#include <stdexcept>

#include "Integral.hpp"

Integral::Integral(int p_, int order_, 
    bool is_free_molecular,
    const SimpleSection* section) : 
        p(p_), order(order_), 
        is_free_molecular(is_free_molecular),
        section(section)
{
}

void Integral::collide(double t, 
                       std::vector<Polygon*>& spacemesh, 
                       const std::vector<int>& mypolys, 
                       Gas& gas)
{
    if (!is_free_molecular) {
        if (order == 1)
        {
            gas.ciGen(2*t, p, section);
            for(size_t i = 0; i < mypolys.size(); i++) {
                gas.ciIter(spacemesh[mypolys[i]]->f().g());
                spacemesh[mypolys[i]]->f().equategf();
            }
        }
        else if (order == 2)
        {
            gas.ciGen(2*t, p/2, section);
            for(size_t i = 0; i < mypolys.size(); i++) {
                gas.ciIter(spacemesh[mypolys[i]]->f().f());
            }

            gas.ciGen(2*t, p/2, section);
            for(size_t i = 0; i < mypolys.size(); i++) {
                gas.ciIter(spacemesh[mypolys[i]]->f().f());
                spacemesh[mypolys[i]]->f().meanf();
            }
        }
        else 
        {
            throw std::invalid_argument("wrong order: " + toStr(order));
        }
    }
}

