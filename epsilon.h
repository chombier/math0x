#ifndef MATH0X_EPSILON_H
#define MATH0X_EPSILON_H

#include <limits>
#include <math0x/types.h>

namespace math0x {

	// machine precision
	template<class U>
	constexpr U epsilon() { return std::numeric_limits<RR>::epsilon(); }
	
	constexpr RR epsilon() { return epsilon<RR>(); }
	
}

#endif
