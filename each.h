#ifndef MATH0X_EACH_H
#define MATH0X_EACH_H

#include <math0x/types.h>

namespace math0x {

	// we need to capture all kinds of refs otherwise we screw up
	// mutable lambdas
	template<class What, class F>
	void each(const What& v, F&& f) {
		for(NN i = 0, n = v.size(); i < n; ++i) f(i);
	};

	// FIXME it's broken lol !
	// note: we want an explicit copy of f
	template<class What, class F>
	void omp_each(const What& v, F f) {
		NN n = v.size();
#pragma omp parallel for firstprivate(f)
		for(NN i = 0; i < n; ++i) f(i);
	};

}


#endif
