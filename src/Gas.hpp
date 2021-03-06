#ifndef _GAS_HPP_
#define _GAS_HPP_

#include "property_tree/property_tree.hpp"

#include "macro.hpp"
#include "maxwell.hpp"
#include "grad13.hpp"
#include "distribution_function.hpp"
#include "axis.hpp"
#include "section.hpp"

class Gas {
    public:
        virtual size_t size() const = 0;

        virtual const std::string macro(const DistributionFunction& f) const = 0;
        virtual double density(const DistributionFunction& f) const = 0;

        virtual const DistributionFunction maxwell(const PropertyTree& tree) const = 0;
        virtual const DistributionFunction grad13(const PropertyTree& tree) const = 0;

        virtual double cut() const = 0;
        virtual double dot(const int i, const V3d n) const = 0;
        virtual const std::vector<V3d>& vel() const = 0; 
        virtual int mirror(const int i, const Axis a) const = 0;
        virtual void equateStreams(DistributionFunction& g,
                             const DistributionFunction& f,
                             const V3d n) const = 0;

        virtual void ciGen(const double time_step, const int p,
                const SimpleSection* section) = 0;


        virtual double ciIter(DistributionFunction& f) const = 0;
};

template <Symmetry symmetry, template<Symmetry> class XiMeshType, typename ColliderType>
class GasTemplate : public Gas {
    public:
        typedef XiMeshType<symmetry> XiMesh;

        GasTemplate(const XiMesh& ximesh_) :
                ximesh(ximesh_), collider(ximesh) {}

        size_t size() const {
            return ximesh.size();
        }

        double density(const DistributionFunction& f) const {
            return ::macro(f, ximesh).n;
        }

        const std::string macro(const DistributionFunction& f) const {
            return toStr(::macro(f, ximesh));
        }

        const DistributionFunction maxwell(const PropertyTree& tree) const {
            return Maxwell(tree, ximesh).func();
        }

        const DistributionFunction grad13(const PropertyTree& tree) const {
            return Grad13(tree, ximesh).func();
        }

        double cut() const {
            return ximesh.cut();
        }

        double dot(const int i, const V3d n) const {
            return ::dot(ximesh[i], n);
        }

        const std::vector<V3d>& vel() const {
            return ximesh.vel();
        }

        int mirror(const int i, const Axis a) const {
            return ximesh.mirror(i, a);
        }

        void equateStreams(DistributionFunction& g,
                     const DistributionFunction& f,
                     const V3d n) const {
            return ::equateStreams(g, f, n, ximesh);
        }

        void ciGen(const double time_step, const int p,
                const SimpleSection* section) {
            collider.gen(time_step, p, ximesh, section);
        }

        double ciIter(DistributionFunction& f) const {
            return collider.iter(f);
        }

    private:
        const XiMesh ximesh;
        ColliderType collider;
};

#endif
