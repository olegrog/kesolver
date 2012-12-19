#ifndef _FACETFACTORY_HPP_
#define _FACETFACTORY_HPP_

#include "singleton/singleton.hpp"
#include "factory/factory.hpp"

#include "PhysicalFacet.hpp"

typedef SingletonHolder< Factory<PhysicalFacet> > FacetFactory;

#define REGISTER_FACET(FacetClass, Name)          \
namespace {                                       \
    PhysicalFacet* create ## FacetClass() {       \
        return new FacetClass;                    \
    }                                             \
    bool is_registered =                          \
        FacetFactory::instance().add(Name, create ## FacetClass); \
}

#endif // _FACETFACTORY_HPP_
