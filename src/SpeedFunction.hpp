#ifndef _SPEED_FUNCTION_
#define _SPEED_FUNCTION_

#include <vector>

#include "v.hpp"
#include "distribution_function.hpp"
#include "auxiliary.hpp"

class SpeedFunction {
	public:
        SpeedFunction() : size_(0) {}

        template <typename F>
		void f(const F& f_);

        size_t size() const { return size_; }

		DistributionFunction& f() { return f1; }
		DistributionFunction& g() { return f2; }

		void equategf();
		void meanf();

		DistributionFunction3& getGradient() { return df; }
		void giveMemoryToGradient();

		DistributionFunction& getFMin() { return fmin; }
		DistributionFunction& getFMax() { return fmax; }
		DistributionFunction& getPhi()  { return phi; }

	private:
		size_t size_;
		DistributionFunction  f1;
		DistributionFunction  f2;
		DistributionFunction3 df;
		DistributionFunction  fmin;
		DistributionFunction  fmax;
		DistributionFunction  phi;
};

template <typename F>
void SpeedFunction::f(const F& f_)
{
    size_ = f_.size();
    copy(f_, f1);
    copy(f_, f2);
}

#endif
