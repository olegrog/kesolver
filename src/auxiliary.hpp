#ifndef _AUXILIARY_HPP_
#define _AUXILIARY_HPP_

#include <string>
#include <sstream>
#include <cmath>

template <typename T> std::string toStr(const T& something) {
	std::ostringstream ss;
	ss << something;
	return ss.str();
}

template <typename T> T strTo(const std::string& str) {
	std::istringstream ss(str);
	T res;
	ss >> res;
	return res;
}

template <typename T>
inline int toInt(const T x) { return static_cast<int>(x); }

template <typename T> T sqr(T x) {
	return x*x;
}

template <typename T> T cube(T x) {
	return x*x*x;
}

template <typename T> T cbrt(T x) {
    return std::pow(x, 1./3.);
}

template <typename F1, typename F2>
void copy(const F1& f1, F2& f2)
{
    f2.resize(f1.size());
    for (size_t i = 0; i < f1.size(); ++i)
        f2[i] = f1[i];
}

const double ln2 = 6.931471806E-1;

#define LABEL std::cout << "LABEL in " << __FILE__ << ':' <<  __LINE__ << std::endl;

#endif // _AUXILIARY_HPP_
