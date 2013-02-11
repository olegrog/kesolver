#include <cassert>

#include "SpeedFunction.hpp"
#include "auxiliary.hpp"

void SpeedFunction::giveMemoryToGradient()
{
	df.resize(size_);
	fmin.resize(size_);
	fmax.resize(size_);
	phi.resize(size_);
}

void SpeedFunction::equategf() {
    copy(f2, f1);
}

void SpeedFunction::equatefg() {
    copy(f1, f2);
}

void SpeedFunction::meanf() {
    for (size_t i = 0; i < size_; ++i) {
        f1[i] = 0.5 * (f1[i] + f2[i]);
        f2[i] = f1[i];
    }
}

